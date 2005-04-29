#include "Interaction_Model_sLepton_EW.H"
#include "MathTools.H"
#include "Message.H"
#include "Run_Parameter.H"
#include <stdio.h>


using namespace AMEGIC;
using namespace ATOOLS;
using namespace std;

Interaction_Model_sLepton_EW::Interaction_Model_sLepton_EW(MODEL::Model_Base * _model,
							   std::string _cplscheme,std::string _yukscheme) :
  Interaction_Model_Base(_model,_cplscheme,_yukscheme)
{ 
  double Ecms2 = sqr(rpa.gen.Ecms());

  g1       = Kabbala(string("g_1"),
		     sqrt(4.*M_PI*ScalarFunction(string("alpha_QED"),Ecms2)));
  g2       = Kabbala(string("g_1/\\sin\\theta_W"), 
		     g1.Value()/sqrt(ScalarConstant(string("sin2_thetaW"))));
  sintW    = Kabbala(string("\\sin\\theta_W"),
		     sqrt(ScalarConstant(string("sin2_thetaW"))));
  costW    = Kabbala(string("\\cos\\theta_W"),
		     sqrt(1.-ScalarConstant(string("sin2_thetaW"))));
  PL       = Kabbala(string("P_L"),1.);
  PR       = Kabbala(string("P_R"),1.);
  M_I      = Kabbala(string("i"),Complex(0.,1.));
  vev      = Kabbala(string("v_{EW}"),ScalarConstant(string("vev")));
  
  v1     = Kabbala(string("v_1"), vev.Value() *
		   sqrt(1./(1.+sqr(ScalarConstant(std::string("tan(beta)"))))));
  v2     = Kabbala(string("v_2"),v1.Value()*ScalarConstant(std::string("tan(beta)")));

  mu       = Kabbala(string("h"),ScalarConstant(string("mu")));
  conj_mu  = Kabbala(string("h"),ScalarConstant(string("mu")));
  K_zero   = Kabbala(string("zero"),0.);
  num_2    = Kabbala(string("2"),2.);    	
  num_3    = Kabbala(string("3"),3.);    	
  num_4    = Kabbala(string("4"),4.);    		
  root2    = Kabbala(string("\\sqrt{2}"),sqrt(2.));
  invroot2 = Kabbala(string("1/\\sqrt{2}"),sqrt(.5));
}

void Interaction_Model_sLepton_EW::c_SSS(Single_Vertex* vertex,int& vanz)
{
  Kabbala kcpl0,kcpl1;

  //sneutrino - Higgs - sneutrino
  for (short int i=81;i<84;i++) {
    Flavour flav = Flavour(kf::code(i));
    for (short int k=31;k<33;k++) {
      Flavour flh = Flavour(kf::code(k));
      if (flh.IsOn() && flav.IsOn()) {
	vertex[vanz].in[0] = flav;
	vertex[vanz].in[1] = flh;
	vertex[vanz].in[2] = flav;

	kcpl0 = -M_I*g2*g2/(costW*costW*num_4)*K_B_R(k-31);
	kcpl1 = kcpl0;
	
	vertex[vanz].cpl[0]  = kcpl0.Value(); 
	vertex[vanz].cpl[1]  = kcpl1.Value();
	vertex[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();
	vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;

	vertex[vanz].ncf   = 1;
	vertex[vanz].Color = new Color_Function(cf::None); 
	
	vertex[vanz].nlf     = 1;
	vertex[vanz].Lorentz = new Lorentz_Function(lf::SSS);
	
	vertex[vanz].on      = 1;
	vanz++;
	//checked RK & FK
      }
    }
  }

  //sneutrino - Hmin - slepton
 
  Flavour flHm = Flavour(kf::Hmin);
  if (flHm.IsOn()) {
    for (short int i=81;i<84;i++) {
      Flavour flav1 = Flavour(kf::code(i));
      for (short int j=71;j<77;j++) {
	Flavour flav2 =Flavour(kf::code(j));
	if(flav1.IsOn() && flav2.IsOn()){
	  vertex[vanz].in[0] = flav1;
	  vertex[vanz].in[1] = flHm.Bar();
	  vertex[vanz].in[2] = flav2;
     
	  Kabbala K_lI = Kabbala(string("\\frac{\\m M_{")+flav2.TexName()+string("}}{ v_1}\\sqrt{2}"),
				 -Flavour(kf::code(2*gen_sLep(flav2)+11)).Yuk()/v1.Value()*sqrt(2.));

	  kcpl0 = M_I*K_Z_Nu(gen_sLep(flav2),i-81)*
	    (-root2*g2*g2/num_4*(v1*K_Z_H(0,0)+v2*K_Z_H(1,0))*
	     K_Z_L(gen_sLep(flav2),j-71)-
	     K_Z_H(0,0)*(K_yuk(flav2)*K_lI*K_Z_L(gen_sLep(flav2),j-71)-
			 (K_l_S(gen_sLep(flav2),0)*K_Z_L(3,j-71)+
			  K_l_S(gen_sLep(flav2),1)*K_Z_L(4,j-71)+
			  K_l_S(gen_sLep(flav2),2)*K_Z_L(5,j-71)))
	     +K_Z_H(1,0)*(K_k_S(gen_sLep(flav2),0)*K_Z_L(3,j-71)+
			  K_k_S(gen_sLep(flav2),1)*K_Z_L(4,j-71)+
			  K_k_S(gen_sLep(flav2),2)*K_Z_L(5,j-71)-
			  K_lI*mu*K_Z_L(gen_sLep(flav2)+3,j-71)));
	 
	  kcpl1 = kcpl0;
	  
	  vertex[vanz].cpl[0]  = kcpl0.Value(); 
	  vertex[vanz].cpl[1]  = kcpl1.Value();
	  vertex[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();
	  vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;
	
	  vertex[vanz].ncf   = 1;
	  vertex[vanz].Color = new Color_Function(cf::None); 
	      
	  vertex[vanz].nlf     = 1;
	  vertex[vanz].Lorentz = new Lorentz_Function(lf::SSS);
	  
	  vertex[vanz].on      = 1;
	  vanz++;
	  
	}
      }
    }
  }
  
  //slepton - A0 - slepton
  Flavour flA0 = Flavour(kf::A0);
  if (flA0.IsOn()) {
    for (short int i=71;i<77;i++) {
      Flavour flav1 = Flavour(kf::code(i));
      for (short int j=71;j<77;j++) {
	Flavour flav2 =Flavour(kf::code(j));
	if(flav1.IsOn() && flav2.IsOn() && i<=j){
	  vertex[vanz].in[0] = flav1;
	  vertex[vanz].in[1] = flA0;
	  vertex[vanz].in[2] = flav2;

	  Kabbala K_lI = Kabbala(string("\\frac{(\\m M_{")+flav1.TexName()+string("})}{ v_1}\\sqrt{2}"),
				 -Flavour(kf::code(2*gen_sLep(flav1)+11)).Yuk()/(v1).Value()*sqrt(2.));
     
	  kcpl0 = (K_l_S(gen_sLep(flav1),gen_sLep(flav2))*
		   (K_Z_L(gen_sLep(flav1),j-71)*K_Z_L(gen_sLep(flav2)+3,i-71)-
		    K_Z_L(gen_sLep(flav1),i-71)*K_Z_L(gen_sLep(flav2)+3,j-71))*K_Z_H(0,0)+
		   K_k_S(gen_sLep(flav1),gen_sLep(flav2))*
		   (K_Z_L(gen_sLep(flav1),j-71)*K_Z_L(gen_sLep(flav2)+3,i-71)-
		    K_Z_L(gen_sLep(flav1),i-71)*K_Z_L(gen_sLep(flav2)+3,j-71))*K_Z_H(1,0)+
		   K_lI*mu*(K_Z_L(gen_sLep(flav1),i-71)*K_Z_L(gen_sLep(flav1)+3,j-71)-
			    K_Z_L(gen_sLep(flav1),j-71)*K_Z_L(gen_sLep(flav1)+3,i-71))*
		   K_Z_H(1,0))*(-invroot2);
	  
	  kcpl1 = kcpl0;
	  
	  vertex[vanz].cpl[0]  = kcpl0.Value(); 
	  vertex[vanz].cpl[1]  = kcpl1.Value();
	  vertex[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();
	  vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;
	  
	  vertex[vanz].ncf   = 1;
	  vertex[vanz].Color = new Color_Function(cf::None); 
	  
	  vertex[vanz].nlf     = 1;
	  vertex[vanz].Lorentz = new Lorentz_Function(lf::SSS);
	  	  
	  vertex[vanz].on      = 1;
	  vanz++;
  
	}
      }
    }
  }
			       
  //slepton - h0/H0 - slepton
 
  for (short int k=31;k<33;k++) {
    Flavour flH = Flavour(kf::code(k)); 
    for (short int i=71;i<77;i++) {
      Flavour flav1 = Flavour(kf::code(i));
      for (short int j=71;j<77;j++) {
	Flavour flav2 =Flavour(kf::code(j));
	if(flH.IsOn() && flav1.IsOn() && flav2.IsOn() && gen_sLep(flav1)==gen_sLep(flav2) && i<=j){
	  vertex[vanz].in[0] = flav1;
	  vertex[vanz].in[1] = flH;
	  vertex[vanz].in[2] = flav2;

	  Kabbala help = K_zero;

	  Kabbala K_lI = Kabbala(string("\\frac{\\m M_{")+flav1.TexName()+string("}}{ v_1}*\\sqrt{2}"),
				 -Flavour(kf::code(2*gen_sLep(flav1)+11)).Yuk()/(v1).Value()*sqrt(2.));
	 
	  Kabbala fac = Kabbala(string("\\frac{1-4sin^2\\theta_W}{2sin^2\\theta_W}"),
				(1.-4.*(sintW).Value()*(sintW).Value())/
				(2.*(sintW).Value()*(sintW).Value()));

	  if (i==j) {help = Kabbala(string("1"),1.);}
 
	  kcpl0 = M_I*(g1*g1/(costW*costW*num_2)*K_B_R(k-31)*
		       (help+fac*K_Z_L(gen_sLep(flav1),i-71)*K_Z_L(gen_sLep(flav1),j-71))
		       - K_lI*K_lI*v1*K_Z_R(0,k-31)*
		       (K_Z_L(gen_sLep(flav1),i-71)*K_Z_L(gen_sLep(flav1),j-71)+
		 	K_Z_L(gen_sLep(flav1)+3,i-71)*K_Z_L(gen_sLep(flav1)+3,j-71))
		     - K_Z_R(0,k-31)/root2*K_l_S(gen_sLep(flav1),gen_sLep(flav2))*
						  (K_Z_L(gen_sLep(flav1),j-71)*
						   K_Z_L(gen_sLep(flav2)+3,i-71)+
						   K_Z_L(gen_sLep(flav1),i-71)*
						   K_Z_L(gen_sLep(flav2)+3,j-71))
		       + K_Z_R(1,k-31)/root2*K_k_S(gen_sLep(flav1),gen_sLep(flav2))*
						  (K_Z_L(gen_sLep(flav1),j-71)*
						   K_Z_L(gen_sLep(flav2)+3,i-71)+
						   K_Z_L(gen_sLep(flav1),i-71)*
						   K_Z_L(gen_sLep(flav2)+3,j-71))
		       - K_lI*K_Z_R(1,k-31)*mu/root2*(K_Z_L(gen_sLep(flav1),i-71)*
						      K_Z_L(gen_sLep(flav1)+3,j-71)+
						      K_Z_L(gen_sLep(flav1),j-71)*
						      K_Z_L(gen_sLep(flav1)+3,i-71)));
						    

	  kcpl1 = kcpl0;

	  vertex[vanz].cpl[0]  = kcpl0.Value();
	  vertex[vanz].cpl[1]  = kcpl1.Value();
	  vertex[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();
	  vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;
	  	  
	  vertex[vanz].ncf   = 1;
	  vertex[vanz].Color = new Color_Function(cf::None); 
	  
	  vertex[vanz].nlf     = 1;
	  vertex[vanz].Lorentz = new Lorentz_Function(lf::SSS);
	  
	  vertex[vanz].on      = 1;
	  vanz++;
  	}
      }
    }
  }
}

void Interaction_Model_sLepton_EW::c_SSV(Single_Vertex* vertex,int& vanz)
{
  Kabbala kcpl0,kcpl1;

  //slepton - Photon - slepton
  Flavour flPh = Flavour(kf::photon);
  if (flPh.IsOn()) {
    for (short int i=71;i<77;i++) {
      Flavour flav = Flavour(kf::code(i));
      Kabbala charge1 = Kabbala(string("Q_{")+flav.TexName()+string("}"),flav.Charge());
      if (flav.IsOn()) {
	vertex[vanz].in[0] = flav;
	vertex[vanz].in[1] = flPh;
	vertex[vanz].in[2] = flav;
	
	kcpl0 = -M_I*g1*charge1;
	kcpl1 = kcpl0;
	
	vertex[vanz].cpl[0]  = kcpl0.Value(); 
	vertex[vanz].cpl[1]  = kcpl1.Value();
	vertex[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();
	vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;

	vertex[vanz].ncf   = 1;
	vertex[vanz].Color = new Color_Function(cf::None); 
	
	vertex[vanz].nlf     = 1;
	vertex[vanz].Lorentz = new Lorentz_Function(lf::SSV);
	vertex[vanz].Lorentz->SetParticleArg(0,2,1);     	

	vertex[vanz].on      = 1;
	vanz++;
      }
    }
  }
 
 //slepton - Z - slepton
  Flavour flZ = Flavour(kf::Z);
  if (flZ.IsOn()) {
    for (short int i=71;i<77;i++) {
      Flavour flav1 = Flavour(kf::code(i));
      for (short int j=i;j<77;j++) {
	Flavour flav2 = Flavour(kf::code(j));
	if (flav1.IsOn() && flav2.IsOn() && gen_sLep(flav1)==gen_sLep(flav2)) {
	  
	  Kabbala help = K_zero;
	  if(i==j) help = sintW*sintW;
	  
	  vertex[vanz].in[0] = flav1;
	  vertex[vanz].in[1] = flZ;
	  vertex[vanz].in[2] = flav2;

	  kcpl0 = M_I*g2/costW*
	    ((K_Z_L(gen_sLep(flav2),j-71)*K_Z_L(gen_sLep(flav2),i-71))/num_2 - help);
	  kcpl1 = kcpl0;	  	  

	  vertex[vanz].cpl[0]  = kcpl0.Value(); 
	  vertex[vanz].cpl[1]  = kcpl1.Value();
	  vertex[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();
	  vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;

	  vertex[vanz].ncf   = 1;
	  vertex[vanz].Color = new Color_Function(cf::None); 
	
	  vertex[vanz].nlf     = 1;
	  vertex[vanz].Lorentz = new Lorentz_Function(lf::SSV);
	  vertex[vanz].Lorentz->SetParticleArg(0,2,1);     	
	  
	  vertex[vanz].on      = 1;
	  vanz++;
	  //checked SS
	}
      }

    }

    //sneutrino - Z - sneutrino
    for (short int i=81;i<84;i++) {
      Flavour flav = Flavour(kf::code(i));
      if (flav.IsOn()) {
 	  vertex[vanz].in[0] = flav;
	  vertex[vanz].in[1] = flZ;
	  vertex[vanz].in[2] = flav;

	  kcpl0 = -M_I*g2/(costW*num_2);
	  kcpl1 = kcpl0;
	  
	  vertex[vanz].cpl[0]  = kcpl0.Value(); 
	  vertex[vanz].cpl[1]  = kcpl1.Value();
	  vertex[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();
	  vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;
	  
	  vertex[vanz].ncf   = 1;
	  vertex[vanz].Color = new Color_Function(cf::None); 
	  
	  vertex[vanz].nlf     = 1;
	  vertex[vanz].Lorentz = new Lorentz_Function(lf::SSV);
	  vertex[vanz].Lorentz->SetParticleArg(0,2,1);     	
	  
	  vertex[vanz].on      = 1;
	  vanz++;
      }
    }
    
  }
  
  //sneutrino - W - slepton
  Flavour flW = Flavour(kf::W);
  if (flW.IsOn()) {
      for (short int i=81;i<84;i++) {
      Flavour flav1 = Flavour(kf::code(i));
      for (short int j=71;j<77;j++) {
	Flavour flav2 =Flavour(kf::code(j));
	if(flav1.IsOn() && flav2.IsOn() && gen_sLep(flav1)==gen_sLep(flav2)) {
	  vertex[vanz].in[0] = flav1;
	  vertex[vanz].in[1] = flW.Bar();
	  vertex[vanz].in[2] = flav2;
	
	  kcpl0 = -M_I*g2/root2
	    *K_Z_Nu(gen_sLep(flav1),gen_sLep(flav1))
	    *K_Z_L(gen_sLep(flav1),j-71);
	  kcpl1 = kcpl0;
	  
	  vertex[vanz].cpl[0]  = kcpl0.Value(); 
	  vertex[vanz].cpl[1]  = kcpl1.Value();
	  vertex[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();
	  vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;
	  
	  vertex[vanz].ncf   = 1;
	  vertex[vanz].Color = new Color_Function(cf::None); 
	  
	  vertex[vanz].nlf     = 1;
	  vertex[vanz].Lorentz = new Lorentz_Function(lf::SSV);
	  vertex[vanz].Lorentz->SetParticleArg(0,2,1);     	
	  
	  vertex[vanz].on      = 1;
	  vanz++;
  	}
      }
    }
  }
}

void Interaction_Model_sLepton_EW::c_SSSS(Single_Vertex* vertex,int& vanz)
{
  Flavour flHmin(kf::Hmin);    
  Flavour flA0(kf::A0);    

  Kabbala kcpl0,kcpl1,help,K_lI, num_1;
  num_1    = Kabbala(string("1"),1.);    	

  for (short int i=71;i<77;++i) {
      Flavour flav1 = Flavour(kf::code(i));
      for (short int j=i;j<77;++j) {
	Flavour flav2 = Flavour(kf::code(j));
	if (flav1.IsOn() && flav2.IsOn() && gen_sLep(flav1)==gen_sLep(flav2)) {

	  //Hmin -> slepton - slepton - Hmin 
	  if (flHmin.IsOn()) {
	    vertex[vanz].in[0] = flHmin;
	    vertex[vanz].in[1] = flav1.Bar();
	    vertex[vanz].in[2] = flav2;
	    vertex[vanz].in[3] = flHmin;
	    
	    vertex[vanz].nleg  = 4;  
	    
	    help = K_zero;
	    if(i==j) help = num_1;
	    
	    K_lI = Kabbala(string("\\frac{\\m M_{")+flav1.TexName()+string("}}{ v_1}*\\sqrt{2}"),
			   -Flavour(kf::code(2*gen_sLep(flav1)+11)).Yuk()/(v1).Value()*sqrt(2.));
	    
	    
	    kcpl0 = M_I*(g1*g1/(costW*costW*num_2)*K_A_H(0,0)*
			 (help-(num_1+num_2*sintW*sintW)/(num_2*sintW*sintW)*
			  K_Z_L(gen_sLep(flav1),i-71)*K_Z_L(gen_sLep(flav1),j-71))-
			 K_lI*K_lI*K_Z_H(0,0)*K_Z_H(0,0)*
			 K_Z_L(gen_sLep(flav1)+3,i-71)*K_Z_L(gen_sLep(flav1)+3,j-71));
	    
	    kcpl1 = kcpl0;
	    
	    vertex[vanz].cpl[0]  = kcpl0.Value(); 
	    vertex[vanz].cpl[1]  = kcpl1.Value();
	    vertex[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();
	    vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;
	    
	    vertex[vanz].ncf   = 1;
	    vertex[vanz].Color = new Color_Function(cf::None); 
	    
	    vertex[vanz].nlf     = 1;
	    vertex[vanz].Lorentz = new Lorentz_Function(lf::SSSS);     
	    
	    vertex[vanz].on      = 1;
	    if (kcpl0.Value()!=Complex(0.,0.) && kcpl1.Value()!=Complex(0.,0.)) vanz++;
	  }
	  //A0 -> slepton - slepton - A0 
	  if (flA0.IsOn()) {
	    vertex[vanz].in[0] = flA0;
	    vertex[vanz].in[1] = flav1.Bar();
	    vertex[vanz].in[2] = flav2;
	    vertex[vanz].in[3] = flA0;
	    
	    vertex[vanz].nleg  = 4;  
	    
	    help = K_zero;
	    if(i==j) help = num_1;
	    
	    K_lI = Kabbala(string("\\frac{\\m M_{")+flav1.TexName()+string("}}{ v_1}*\\sqrt{2}"),
			   -Flavour(kf::code(2*gen_sLep(flav1)+11)).Yuk()/(v1).Value()*sqrt(2.));
	    
	    
	    kcpl0 = M_I*(g1*g1/(costW*costW*num_2)*K_A_H(0,0)*
			 (help+(num_1-num_4*sintW*sintW)/(num_2*sintW*sintW)*
			  K_Z_L(gen_sLep(flav1),i-71)*K_Z_L(gen_sLep(flav1),j-71))-
			 K_lI*K_lI*K_Z_H(0,0)*K_Z_H(0,0)*
			 (K_Z_L(gen_sLep(flav1),i-71)*K_Z_L(gen_sLep(flav1),j-71) +
			 K_Z_L(gen_sLep(flav1)+3,i-71)*K_Z_L(gen_sLep(flav1)+3,j-71)));
	    
	    kcpl1 = kcpl0;
	    
	    vertex[vanz].cpl[0]  = kcpl0.Value(); 
	    vertex[vanz].cpl[1]  = kcpl1.Value();
	    vertex[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();
	    vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;
	    
	    vertex[vanz].ncf   = 1;
	    vertex[vanz].Color = new Color_Function(cf::None); 
	    
	    vertex[vanz].nlf     = 1;
	    vertex[vanz].Lorentz = new Lorentz_Function(lf::SSSS);     
	    
	    vertex[vanz].on      = 1;
	    if (kcpl0.Value()!=Complex(0.,0.) && kcpl1.Value()!=Complex(0.,0.)) vanz++;
	  }
	  //h0/H0 -> slepton - slepton - h0/H0 
	  for (int k=31;k<33;k++) {
	    Flavour flh1 = Flavour(kf::code(k));
	    for (int l=k;l<33;l++) {
	      Flavour flh2 = Flavour(kf::code(l));
	      if (flh1.IsOn() && flh2.IsOn()) {
		
		vertex[vanz].in[0] = flh1;
		vertex[vanz].in[1] = flav1.Bar();
		vertex[vanz].in[2] = flav2;
		vertex[vanz].in[3] = flh2;
		
		vertex[vanz].nleg  = 4;  
		
		help = K_zero;
		if(i==j) help = num_1;
		
		K_lI = Kabbala(string("\\frac{\\m M_{")+flav1.TexName()+string("}}{ v_1}*\\sqrt{2}"),
			       -Flavour(kf::code(2*gen_sLep(flav1)+11)).Yuk()/(v1).Value()*sqrt(2.));
		
		
		kcpl0 = M_I*(g1*g1/(costW*costW*num_2)*K_A_R(k-31,l-31)*
			     (help+(num_1-num_4*sintW*sintW)/(num_2*sintW*sintW)*
			      K_Z_L(gen_sLep(flav1),i-71)*K_Z_L(gen_sLep(flav1),j-71))-
			     K_lI*K_lI*K_Z_R(0,k-31)*K_Z_R(0,l-31)*
			     (K_Z_L(gen_sLep(flav1),i-71)*K_Z_L(gen_sLep(flav1),j-71) +
			      K_Z_L(gen_sLep(flav1)+3,i-71)*K_Z_L(gen_sLep(flav1)+3,j-71)));
		
		kcpl1 = kcpl0;
		
		vertex[vanz].cpl[0]  = kcpl0.Value(); 
		vertex[vanz].cpl[1]  = kcpl1.Value();
		vertex[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();
		vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;
		
		vertex[vanz].ncf   = 1;
		vertex[vanz].Color = new Color_Function(cf::None); 
		
		vertex[vanz].nlf     = 1;
		vertex[vanz].Lorentz = new Lorentz_Function(lf::SSSS);     
		
		vertex[vanz].on      = 1;
		if (kcpl0.Value()!=Complex(0.,0.) && kcpl1.Value()!=Complex(0.,0.)) vanz++;
	      }
	    }
	  }
	}
      }
  }

  for (short int i=81;i<84;++i) {
    Flavour flsnu = Flavour(kf::code(i));
    if (flsnu.IsOn()) {
      //h0/H0 -> snu - snub - h0/H0
        for (int k=31;k<33;k++) {
	    Flavour flh1 = Flavour(kf::code(k));
	    for (int l=k;l<33;l++) {
	      Flavour flh2 = Flavour(kf::code(l));
	      if (flh1.IsOn() && flh2.IsOn()) {
		
		vertex[vanz].in[0] = flh1;
		vertex[vanz].in[1] = flsnu.Bar();
		vertex[vanz].in[2] = flsnu;
		vertex[vanz].in[3] = flh2;
		
		vertex[vanz].nleg  = 4;  
		
		kcpl0 = -M_I*g2*g2/(costW*costW*num_4)*K_A_R(k-31,l-31);
		
		kcpl1 = kcpl0;
		
		vertex[vanz].cpl[0]  = kcpl0.Value(); 
		vertex[vanz].cpl[1]  = kcpl1.Value();
		vertex[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();
		vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;
		
		vertex[vanz].ncf   = 1;
		vertex[vanz].Color = new Color_Function(cf::None); 
		
		vertex[vanz].nlf     = 1;
		vertex[vanz].Lorentz = new Lorentz_Function(lf::SSSS);     
		
		vertex[vanz].on      = 1;
		vanz++;
	      }
	    }
	}
        if (flA0.IsOn()) {
	  // A0 -> snu - snub - A0
	  vertex[vanz].in[0] = flA0;
	  vertex[vanz].in[1] = flsnu.Bar();
	  vertex[vanz].in[2] = flsnu;
	  vertex[vanz].in[3] = flA0;
	  
	  vertex[vanz].nleg  = 4;  
	  
	  kcpl0 = -M_I*g2*g2/(costW*costW*num_4)*K_A_H(0,0);
	  
	  kcpl1 = kcpl0;
	  
	  vertex[vanz].cpl[0]  = kcpl0.Value(); 
	  vertex[vanz].cpl[1]  = kcpl1.Value();
	  vertex[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();
	  vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;
	  
	  vertex[vanz].ncf   = 1;
	  vertex[vanz].Color = new Color_Function(cf::None); 
	  
	  vertex[vanz].nlf     = 1;
	  vertex[vanz].Lorentz = new Lorentz_Function(lf::SSSS);     
	  
	  vertex[vanz].on      = 1;
	  vanz++;
	} 
	if (flHmin.IsOn()) {
	  // H- -> snu - snub - H-
	  vertex[vanz].in[0] = flHmin;
	  vertex[vanz].in[1] = flsnu.Bar();
	  vertex[vanz].in[2] = flsnu;
	  vertex[vanz].in[3] = flHmin;
	  
	  Kabbala cot2TW = Kabbala(string("cot2\\theta_W"),
				   (costW.Value()*costW.Value()-sintW.Value()*sintW.Value())/
				   (2.*sintW.Value()*costW.Value()));
	  
	  
	  Flavour lepton = Flavour(kf::code(11+(i-81)*2));
	  
	  Kabbala K_lI = Kabbala(string("\\frac{\\m M_{")+lepton.TexName()+string("}}{ v_1}\\sqrt{2}"),
				 -lepton.Yuk()/v1.Value()*sqrt(2.));
	  
	  
	  vertex[vanz].nleg  = 4;  
	  
	  kcpl0 = M_I*K_Z_Nu(i-81,i-81)*K_Z_Nu(i-81,i-81)*
	    (g1*g1*cot2TW/(num_2*sintW*costW)*K_A_H(0,0) - K_lI*K_lI*K_Z_H(0,0)*K_Z_H(0,0));
	  
	  kcpl1 = kcpl0;
	  
	  vertex[vanz].cpl[0]  = kcpl0.Value(); 
	  vertex[vanz].cpl[1]  = kcpl1.Value();
	  vertex[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();
	  vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;
	  
	  vertex[vanz].ncf   = 1;
	  vertex[vanz].Color = new Color_Function(cf::None); 
	  
	  vertex[vanz].nlf     = 1;
	  vertex[vanz].Lorentz = new Lorentz_Function(lf::SSSS);     
	  
	  vertex[vanz].on      = 1;
	  vanz++;
	} 
	if (flHmin.IsOn()) {
	  for (int j=71;j<77;j++) {
	    Flavour flslep = Flavour(kf::code(j));
	    if (gen_sLep(flsnu)==gen_sLep(flslep)) {
	      // H- -> sLepton - snub - h0/H0
	      for (int k=31;k<33;k++) {
		Flavour flh = Flavour(kf::code(k));
		if (flh.IsOn()) {
		  vertex[vanz].in[0] = flHmin;
		  vertex[vanz].in[1] = flsnu.Bar();
		  vertex[vanz].in[2] = flslep;
		  vertex[vanz].in[3] = flh;
		  
		  Flavour lepton = Flavour(kf::code(11+gen_sLep(flslep)*2));
		  
		  Kabbala K_lI = Kabbala(string("\\frac{\\m M_{")+lepton.TexName()+string("}}{ v_1}\\sqrt{2}"),
					 -lepton.Yuk()/v1.Value()*sqrt(2.));
		  
		  vertex[vanz].nleg  = 4;  
		  
		  kcpl0 = M_I/root2*K_Z_L(gen_sLep(flslep),j-71)*K_Z_Nu(gen_sLep(flslep),gen_sLep(flsnu))*
		    (-g2*g2/num_2*(K_Z_H(0,0)*K_Z_R(0,k-31) + K_Z_H(1,0)*K_Z_R(1,k-31)) 
		     + K_lI*K_lI*K_Z_H(0,0)*K_Z_R(0,k-31));
		  
		  kcpl1 = kcpl0;
	      
		  vertex[vanz].cpl[0]  = kcpl0.Value(); 
		  vertex[vanz].cpl[1]  = kcpl1.Value();
		  vertex[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();
		  vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;
		  
		  vertex[vanz].ncf   = 1;
		  vertex[vanz].Color = new Color_Function(cf::None); 
		  
		  vertex[vanz].nlf     = 1;
		  vertex[vanz].Lorentz = new Lorentz_Function(lf::SSSS);     
		  
		  vertex[vanz].on      = 1;
		  if (kcpl0.Value()!=Complex(0.,0.) && kcpl1.Value()!=Complex(0.,0.)) vanz++;
		}
	      }
	      // H- -> sLepton - snub - A0
	      if (flA0.IsOn()) {
		vertex[vanz].in[0] = flHmin;
		vertex[vanz].in[1] = flsnu.Bar();
		vertex[vanz].in[2] = flslep;
		vertex[vanz].in[3] = flA0;
		
		Flavour lepton = Flavour(kf::code(11+gen_sLep(flslep)*2));
		
		Kabbala K_lI = Kabbala(string("\\frac{\\m M_{")+lepton.TexName()+string("}}{ v_1}\\sqrt{2}"),
				       -lepton.Yuk()/v1.Value()*sqrt(2.));
		
		vertex[vanz].nleg  = 4;  
		
		kcpl0 = -num_1/root2*K_Z_L(gen_sLep(flslep),j-71)*K_Z_Nu(gen_sLep(flslep),gen_sLep(flsnu))*
		  (g2*g2/num_2*K_A_H(0,0) - K_lI*K_lI*K_Z_H(0,0)*K_Z_H(0,0));
		
		kcpl1 = kcpl0;
		
		vertex[vanz].cpl[0]  = kcpl0.Value(); 
		vertex[vanz].cpl[1]  = kcpl1.Value();
		vertex[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();
		vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;
		
		vertex[vanz].ncf   = 1;
		vertex[vanz].Color = new Color_Function(cf::None); 
		
		vertex[vanz].nlf     = 1;
		vertex[vanz].Lorentz = new Lorentz_Function(lf::SSSS);     
		
		vertex[vanz].on      = 1;
		if (kcpl0.Value()!=Complex(0.,0.) && kcpl1.Value()!=Complex(0.,0.)) vanz++;
	      }
	    }
	  }
	}
    }
  }
  //snu - snu - snu - snu
  for (short int i=81;i<84;++i) {
    Flavour snu1 = Flavour(kf::code(i));
    if (snu1.IsOn()) {
      for (short int j=i;j<84;++j) {
	Flavour snu2 = Flavour(kf::code(j));
	if (snu2.IsOn()) {
	  
	  vertex[vanz].in[0] = snu1;
	  vertex[vanz].in[1] = snu1;
	  vertex[vanz].in[2] = snu2;
	  vertex[vanz].in[3] = snu2.Bar();
	  
	  vertex[vanz].nleg  = 4;  
	  
	  Kabbala fac = num_1;
	  
	  if (i==j) fac = num_2;
	  
	  kcpl0 = -M_I*g2*g2/(num_4*costW*costW)*fac;
      
	  kcpl1 = kcpl0;
      	  
	  vertex[vanz].cpl[0]  = kcpl0.Value(); 
	  vertex[vanz].cpl[1]  = kcpl1.Value();
	  vertex[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();
	  vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;
	  
	  vertex[vanz].ncf   = 1;
	  vertex[vanz].Color = new Color_Function(cf::None); 
	  
	  vertex[vanz].nlf     = 1;
	  vertex[vanz].Lorentz = new Lorentz_Function(lf::SSSS);     
	  
	  vertex[vanz].on      = 1;
	  if (kcpl0.Value()!=Complex(0.,0.) && kcpl1.Value()!=Complex(0.,0.)) vanz++;
	
	  //snu - sLep - snu - sLep

	  for (int k=71;k<77;k++) {
	    Flavour slep1 = Flavour(kf::code(k));
	    if (slep1.IsOn()) {
	      for (int l=k;l<77;l++) {
		Flavour slep2 = Flavour(kf::code(l));
		if (slep2.IsOn()) {
		  
		  vertex[vanz].in[0] = snu1;
		  vertex[vanz].in[1] = snu2;
		  vertex[vanz].in[2] = slep1;
		  vertex[vanz].in[3] = slep2.Bar();
		  
		  vertex[vanz].nleg  = 4;  
		  
		  Kabbala help1 = K_zero;
		  if (k==l) help1 = num_1;

		  Kabbala help2 = K_zero;
		  if (i==j) help2 = num_1;
		  
		  Kabbala addendum = K_zero;

		  for (int m=0;m<3;m++) {
		    for (int n=0;n<3;n++) {
		      addendum += (g2*g2/num_2*K_Z_L(m,l-71)*K_Z_L(n,k-71) + 
				   K_l(m)*K_l(n)*K_Z_L(m+3,l-71)*K_Z_L(n+3,k-71))
			*K_Z_Nu(n,gen_sLep(snu1))*K_Z_Nu(m,gen_sLep(snu2));
		    }
		  }
		  
		  kcpl0 = -M_I*(-g1*g1/(num_2*costW*costW)*
				(help1 + (num_1-num_4*sintW*sintW)/(num_2*sintW*sintW)*
				 (K_Z_L(0,l-71)*K_Z_L(0,k-71)+
				  K_Z_L(1,l-71)*K_Z_L(1,k-71)+
				  K_Z_L(2,l-71)*K_Z_L(2,k-71)))*help2 + addendum);
				
		  kcpl1 = kcpl0;
		  
		  vertex[vanz].cpl[0]  = kcpl0.Value(); 
		  vertex[vanz].cpl[1]  = kcpl1.Value();
		  vertex[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();
		  vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;
		  
		  vertex[vanz].ncf   = 1;
		  vertex[vanz].Color = new Color_Function(cf::None); 
		  
		  vertex[vanz].nlf     = 1;
		  vertex[vanz].Lorentz = new Lorentz_Function(lf::SSSS);     
		  
		  if (kcpl0.Value()!=Complex(0.,0.) && kcpl1.Value()!=Complex(0.,0.)) vertex[vanz].on      = 1;
		  vanz++;
		}
	      }
	    }
	  }
	}
      }
    }
  }
  //slep - slep - slep - slep 
  for (int i=71;i<77;i++) {
    Flavour slepi = Flavour(kf::code(i));
    if (slepi.IsOn()) {
      for (int j=71;j<77;j++) {
	Flavour slepj = Flavour(kf::code(j));
	if (slepj.IsOn()) {
	  for (int k=71;k<77;k++) {
	    Flavour slepk = Flavour(kf::code(k));
	    if (slepk.IsOn()) {
	      for (int l=71;l<77;l++) {
		Flavour slepl = Flavour(kf::code(l));
		if (slepl.IsOn()) {
	
		  vertex[vanz].in[0] = slepi;
		  vertex[vanz].in[1] = slepk;
		  vertex[vanz].in[2] = slepl;
		  vertex[vanz].in[3] = slepj.Bar();
		  
		  vertex[vanz].nleg  = 4;  
		  
		  Kabbala d_il = K_zero;
		  if (i==l) d_il = num_1;
		  
		  Kabbala d_jk = K_zero;
		  if (j==k) d_jk = num_1;
		  
		  Kabbala d_ik = K_zero;
		  if (i==k) d_ik = num_1;
		  
		  Kabbala d_jl = K_zero;
		  if (j==l) d_jl = num_1;
		  
		  Kabbala ZLsum_il = K_zero;
		  Kabbala ZLsum_ik = K_zero;
		  Kabbala ZLsum_jl = K_zero;
		  Kabbala ZLsum_jk = K_zero;
		 
		  for (int t=0;t<3;t++) {
		    ZLsum_il += K_Z_L(t,i-71)*K_Z_L(t,l-71);
		    ZLsum_ik += K_Z_L(t,i-71)*K_Z_L(t,k-71);
		    ZLsum_jl += K_Z_L(t,j-71)*K_Z_L(t,l-71);
		    ZLsum_jk += K_Z_L(t,j-71)*K_Z_L(t,k-71);
		  }
		  
		  Kabbala addendum = K_zero;
		  
		  for (int m=0;m<3;m++) {
		    for (int n=0;n<3;n++) {
		      addendum += K_l(m)*K_l(n)*
			(K_Z_L(m+3,i-71)*K_Z_L(n,j-71) + K_Z_L(m+3,j-71)*K_Z_L(n,i-71))*
			(K_Z_L(m,l-71)*K_Z_L(n+3,k-71) + K_Z_L(m,k-71)*K_Z_L(n+3,l-71)); 
		    }
		  }
		  
		  kcpl0 = -M_I*(g1*g1/(costW*costW)*(d_il*d_jk+d_ik*d_jl)
				- num_3*g1*g1/(num_2*costW*costW)*
				(d_il*ZLsum_jk+d_ik*ZLsum_jl+d_jl*ZLsum_ik+d_jk*ZLsum_il) +
				g2*g2*(num_1+num_2*num_4*sintW*sintW)/(num_4*costW*costW)*
				(ZLsum_il*ZLsum_jl + ZLsum_ik*ZLsum_jl) + 
				addendum);
		  
		  kcpl1 = kcpl0;
		  
		  vertex[vanz].cpl[0]  = kcpl0.Value(); 
		  vertex[vanz].cpl[1]  = kcpl1.Value();
		  vertex[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();
		  vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;
		  
		  vertex[vanz].ncf   = 1;
		  vertex[vanz].Color = new Color_Function(cf::None); 
		  
		  vertex[vanz].nlf     = 1;
		  vertex[vanz].Lorentz = new Lorentz_Function(lf::SSSS);     
		  
		  vertex[vanz].on      = 1;
		  if (kcpl0.Value()!=Complex(0.,0.) && kcpl1.Value()!=Complex(0.,0.)) vanz++;
		}
	      }
	    }
	  }
	}
      }
    }
  }
}


void Interaction_Model_sLepton_EW::c_SSVV(Single_Vertex* vertex,int& vanz)
{
  Flavour flavW(kf::W);
  Flavour flavZ(kf::Z);
  Flavour flavPhoton(kf::photon);
  Kabbala kcpl0,kcpl1,help,num_1;
  num_1    = Kabbala(string("1"),1.);    	


  for (short int i=71;i<77;i++) {
    Flavour flav = Flavour(kf::code(i));
    if(flav.IsOn() && flavPhoton.IsOn()){
      
      // P - L - L - P  
      vertex[vanz].in[0] = flavPhoton;
      vertex[vanz].in[1] = flav.Bar();
      vertex[vanz].in[2] = flav;
      vertex[vanz].in[3] = flavPhoton;
      
      vertex[vanz].nleg     = 4;
      
      kcpl0 = num_2*M_I*g1*g1;
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
  for (short int i=71;i<77;i++) {
    Flavour flav1 = Flavour(kf::code(i));
    for (short int j=71;j<77;j++) {
      Flavour flav2 = Flavour(kf::code(j));
      if(flav1.IsOn() && flav2.IsOn() && gen_sLep(flav1)==gen_sLep(flav2)) {
	if(flavZ.IsOn()) {
	  help = K_zero;
	  if(i==j) help = num_1;
	  if(flavPhoton.IsOn()) {
	    // P - L - L - Z
	    vertex[vanz].in[0] = flavPhoton;
	    vertex[vanz].in[1] = flav1.Bar();
	    vertex[vanz].in[2] = flav2;
	    vertex[vanz].in[3] = flavZ;
	    
	    vertex[vanz].nleg     = 4;
	    
	    kcpl0 = M_I*g1*g2/costW*((K_Z_L(0,i-71)*K_Z_L(0,j-71) + 
				      K_Z_L(1,i-71)*K_Z_L(1,j-71) +
				      K_Z_L(2,i-71)*K_Z_L(2,j-71))
				      - num_2*help*sintW*sintW);
	    
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
	  // Z - L - L - Z
	    vertex[vanz].in[0] = flavZ;
	    vertex[vanz].in[1] = flav1.Bar();
	    vertex[vanz].in[2] = flav2;
	    vertex[vanz].in[3] = flavZ;
	    
	    vertex[vanz].nleg     = 4;
	    
	    kcpl0 = num_2*M_I*g1*g1/(costW*costW)*(help*sintW*sintW +
						   (num_1 - num_4*sintW*sintW)/(num_4*sintW*sintW)*
						   (K_Z_L(0,i-71)*K_Z_L(0,j-71) +
						    K_Z_L(1,i-71)*K_Z_L(1,j-71) +
						    K_Z_L(2,i-71)*K_Z_L(2,j-71)));
	    
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
  }
  
  for (short int i=81;i<84;i++) {
    Flavour flav = Flavour(kf::code(i));
    if(flav.IsOn()) { 
      
      if (flavZ.IsOn()) {
	// Z - snu - snu - Z  
	vertex[vanz].in[0] = flavZ;
	vertex[vanz].in[1] = flav.Bar();
	vertex[vanz].in[2] = flav;
	vertex[vanz].in[3] = flavZ;
	
	vertex[vanz].nleg     = 4;
	
	kcpl0 = M_I*g2*g2/(costW*costW*num_2);
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
      if (flavW.IsOn()) {
	// W - snu - snu - W 
	vertex[vanz].in[0] = flavW;
	vertex[vanz].in[1] = flav.Bar();
	vertex[vanz].in[2] = flav;
	vertex[vanz].in[3] = flavW;
	
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
  for (short int i=71;i<77;i++) {
    Flavour flav1 = Flavour(kf::code(i));
    for (short int j=i;j<77;j++) {
      Flavour flav2 =Flavour(kf::code(j));
      if(flav1.IsOn() && flav2.IsOn() && gen_sLep(flav1)==gen_sLep(flav2)){
	
	// W - L - L - W  
	if (flavW.IsOn() && gen_sLep(flav1)==gen_sLep(flav2)) {
	  
	  vertex[vanz].in[0] = flavW;
	  vertex[vanz].in[1] = flav1.Bar();
	  vertex[vanz].in[2] = flav2;
	  vertex[vanz].in[3] = flavW;
	  
	  vertex[vanz].nleg     = 4;
	  
	  kcpl0 = M_I*g2*g2/num_2*K_Z_L(gen_sLep(flav1),i-71)*K_Z_L(gen_sLep(flav1),j-71);
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
  }
  
  // W - snu - L - P/Z  
  for (short int i=81;i<84;i++) {
    Flavour flav1 = Flavour(kf::code(i));
    for (short int j=71;j<77;j++) {
      Flavour flav2 =Flavour(kf::code(j));
      if(flav1.IsOn() && flav2.IsOn() && gen_sLep(flav1)==gen_sLep(flav2)){
	
	// W - snu - L - P  
	if (flavW.IsOn()) {
	  if (flavPhoton.IsOn()) {
	    
	    vertex[vanz].in[0] = flavW;
	    vertex[vanz].in[1] = flav1.Bar();
	    vertex[vanz].in[2] = flav2;
	    vertex[vanz].in[3] = flavPhoton;
	    
	    vertex[vanz].nleg     = 4;
	    
	    kcpl0 = -M_I*g1*g2/root2*(K_Z_Nu(gen_sLep(flav1),i-81)*
				      K_Z_L(gen_sLep(flav1),j-71));
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
	  
	  if (flavZ.IsOn()) {
	    // W - snu - L - Z  
	    
	    vertex[vanz].in[0] = flavW;
	    vertex[vanz].in[1] = flav1.Bar();
	    vertex[vanz].in[2] = flav2;
	    vertex[vanz].in[3] = flavZ;
	    
	    vertex[vanz].nleg     = 4;
	    
	    kcpl0 = M_I*g1*g1/(root2*costW)*(K_Z_Nu(gen_sLep(flav1),i-81)*
					     K_Z_L(gen_sLep(flav1),j-71));
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
    }
  }
}

void Interaction_Model_sLepton_EW::c_FFS(Single_Vertex* vertex,int& vanz){

  Kabbala kcpl0,kcpl1;
  
  Kabbala K_l1 = Kabbala(string("l^1"),
			 -Flavour(kf::code(11)).Yuk()/v1.Value()*sqrt(2.));

  Kabbala K_l2 = Kabbala(string("l^2"),
			 -Flavour(kf::code(13)).Yuk()/v1.Value()*sqrt(2.));
  
  Kabbala K_l3 = Kabbala(string("l^3"),
			 -Flavour(kf::code(15)).Yuk()/v1.Value()*sqrt(2.));
  
  //neutrino - sneutrino - neutralino
  
  for (short int i=12;i<17;i+=2) {
    Flavour flav1 = Flavour(kf::code(i));
    for (short int j=43;j<47;j++) {
      Flavour flav2 = Flavour(kf::code(j));
      for (short int k=81;k<84;k++) {
	Flavour flav3 = Flavour(kf::code(k));
	if (flav1.IsOn() && flav2.IsOn() && flav3.IsOn()) {
	  vertex[vanz].in[0] = flav1;
	  vertex[vanz].in[1] = flav3;
	  vertex[vanz].in[2] = flav2;
	  
	  kcpl0 = K_zero;
	  kcpl1 = M_I*g2/(costW*root2)*K_Z_Nu((i-12)/2,k-81)*
	    (K_Z_N(0,j-43)*sintW-K_Z_N(1,j-43)*costW);
	  
	  vertex[vanz].cpl[0] = kcpl0.Value();
	  vertex[vanz].cpl[1] = kcpl1.Value();
	  vertex[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();
	  vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;
	
	  vertex[vanz].ncf   = 1;
	  vertex[vanz].Color = new Color_Function(cf::None); 
	  
	  vertex[vanz].nlf     = 1;
	  vertex[vanz].Lorentz = new Lorentz_Function(lf::FFS);
	  
	  vertex[vanz].on      = 1;
	  vanz++;
	  //checked SS
	}
      }
    }
  }  
    
  //neutrino - slepton   - Chargino  

  for (short int i=12;i<17;i+=2) {
    Flavour flav1 = Flavour(kf::code(i));
    
    Kabbala K_lI = Kabbala(string("\\frac{\\m M_{")+flav1.TexName()+string("}}{ v_1}\\sqrt{2}"),
			   -Flavour(kf::code(i-1)).Yuk()/v1.Value()*sqrt(2.));
   
    for (short int j=41;j<43;j++) {
      Flavour flav2 = Flavour(kf::code(j));
      for (short int k=71;k<77;k++) {
	Flavour flav3 = Flavour(kf::code(k));
	if (flav1.IsOn() && flav2.IsOn() && flav3.IsOn()) {
	  vertex[vanz].in[0] = flav1;
	  vertex[vanz].in[1] = flav3;
	  vertex[vanz].in[2] = flav2;
	 
	  kcpl0 = K_zero;
	  kcpl1 = K_yuk_sign(flav2)*(-M_I)*(K_Z_MI(0,j-41)*K_Z_L((i-12)/2,k-71)*g2+
					    K_Z_MI(1,j-41)*K_Z_L((i-12)/2+3,k-71)*K_lI);
	  
	  
	  vertex[vanz].cpl[0] = kcpl0.Value();
	  vertex[vanz].cpl[1] = kcpl1.Value();
	  vertex[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();
	  vertex[vanz].cpl[2] = 0.;vertex[vanz].cpl[3]  = 0.;
	  
	  vertex[vanz].ncf   = 1;
	  vertex[vanz].Color = new Color_Function(cf::None); 
	  
	  vertex[vanz].nlf     = 1;
	  vertex[vanz].Lorentz = new Lorentz_Function(lf::FFS);
	  
	  vertex[vanz].on     = 1;
	  vanz++;
	  //checked SS
	}
      }
    }
  }

  //lepton - Chargino - sneutrino  
 
  for (short int i=11;i<16;i+=2) {
    Flavour flav1 = Flavour(kf::code(i));
  
    Kabbala K_lI = Kabbala(string("\\frac{(\\m M_{")+flav1.TexName()+string(")}}{ v_1}\\sqrt{2}"),
			   -flav1.Yuk()/v1.Value()*sqrt(2.));
       
    for (short int j=41;j<43;j++) {
      Flavour flav2 = Flavour(kf::code(j));
      for (short int k=81;k<84;k++) {
	Flavour flav3 = Flavour(kf::code(k));
	if (flav1.IsOn() && flav2.IsOn() && flav3.IsOn()) {
	  vertex[vanz].in[0] = flav1;
	  vertex[vanz].in[1] = flav3;
	  vertex[vanz].in[2] = flav2;

	  kcpl0 = -M_I*K_lI*K_Z_MI(1,j-41)*K_Z_Nu((i-11)/2,k-81);
	  kcpl1 = -M_I*g2*K_Z_PL(0,j-41)*K_Z_Nu((i-11)/2,k-81); 

	  vertex[vanz].cpl[0] = kcpl0.Value();
	  vertex[vanz].cpl[1] = kcpl1.Value();
	  vertex[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();
	  vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;

	  vertex[vanz].ncf   = 1;
	  vertex[vanz].Color = new Color_Function(cf::None); 
	  
	  vertex[vanz].nlf     = 1;
	  vertex[vanz].Lorentz = new Lorentz_Function(lf::FFS);
	  
	  vertex[vanz].on      = 1;
	  vanz++;
	  //checked RK & FK & SS
	}
      }
    }
  }

  //lepton - slepton - Neutralino  
 
  for (short int i=11;i<16;i+=2) {
    Flavour flav1 = Flavour(kf::code(i));
    for (short int j=43;j<47;j++) {
      Flavour flav2 = Flavour(kf::code(j));
      for (short int k=71;k<77;k++) {
	Flavour flav3 = Flavour(kf::code(k));
	if (flav1.IsOn() && flav2.IsOn() && flav3.IsOn()) {
	  
	  //fixed + save
	  vertex[vanz].in[0] = flav1;
	  vertex[vanz].in[1] = flav3;
	  vertex[vanz].in[2] = flav2;

	  Kabbala K_lI = Kabbala(string("\\frac{\\m M_{")+flav1.TexName()+string("}}{ v_1}\\sqrt{2}"),
				     -flav1.Yuk()/v1.Value()*sqrt(2.));
	  
	  kcpl0 = -M_I*(g1/costW*root2*K_Z_L((i-11)/2+3,k-71)*K_Z_N_com_conj(0,j-43)
			-K_lI*K_Z_L((i-11)/2,k-71)*K_Z_N_com_conj(2,j-43));
	  
	  kcpl1 = M_I*(g2/(costW*root2)*K_Z_L((i-11)/2,k-71)*
	    (K_Z_N_com(0,j-43)*sintW + K_Z_N_com(1,j-43)*costW) + K_lI*K_Z_L((i-11)/2+3,k-71)*K_Z_N_com(2,j-43));
	  
	  vertex[vanz].cpl[0] = kcpl0.Value();
	  vertex[vanz].cpl[1] = kcpl1.Value();
	  vertex[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();
	  vertex[vanz].cpl[2] = 0.;vertex[vanz].cpl[3]  = 0.;
	  
	  vertex[vanz].ncf   = 1;
	  vertex[vanz].Color = new Color_Function(cf::None); 
	  
	  vertex[vanz].nlf     = 1;
	  vertex[vanz].Lorentz = new Lorentz_Function(lf::FFS);
	  
	  vertex[vanz].on     = 1;
	  vanz++;
	  //checked RK & SS
	}
      }
    }
  }
}

Kabbala Interaction_Model_sLepton_EW::K_l(short int i) 
{
  char hi[2];
  sprintf(hi,"%i",i);
  
  return Kabbala(string("l^")+string(hi),
		 -Flavour(kf::code(2*i+11)).Yuk()/v1.Value()*sqrt(2.));

}

Kabbala Interaction_Model_sLepton_EW::K_Z_Nu(short int i,short int j)       
{   
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("Z^{")+string(hi)+string(hj)+string("}_\\nu"),
		 ComplexMatrixElement(string("Z_nu"),i,j));
}  
 
Kabbala Interaction_Model_sLepton_EW::K_Z_L(short int i,short int j)       
{   
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("Z^{")+string(hi)+string(hj)+string("}_U"),
		 ComplexMatrixElement(string("Z_l"),i,j));
}  

Kabbala Interaction_Model_sLepton_EW::K_A_H(short int i,short int j) {
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("A^{")+string(hi)+string(hj)+string("}_H"),
		 ComplexMatrixElement(std::string("Z_H"),0,i) *
		 ComplexMatrixElement(std::string("Z_H"),0,j) -
		 ComplexMatrixElement(std::string("Z_H"),1,i) *
		 ComplexMatrixElement(std::string("Z_H"),1,j));
}  

Kabbala Interaction_Model_sLepton_EW::K_A_R(short int i,short int j) {
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("A^{")+string(hi)+string(hj)+string("}_R"),
		 ComplexMatrixElement(std::string("Z_R"),0,i) *
		 ComplexMatrixElement(std::string("Z_R"),0,j) -
		 ComplexMatrixElement(std::string("Z_R"),1,i) *
		 ComplexMatrixElement(std::string("Z_R"),1,j));
}  

Kabbala Interaction_Model_sLepton_EW::K_k_S(short int i,short int j)
{   
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("e^{")+string(hi)+string(hj)+string("}_S"),
		 ComplexMatrixElement(string("ks"),i,j));
}  

Kabbala Interaction_Model_sLepton_EW::K_l_S(short int i,short int j)
{   
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("d^{")+string(hi)+string(hj)+string("}_S"),
		 ComplexMatrixElement(string("ls"),i,j));
}

Kabbala Interaction_Model_sLepton_EW::K_yuk(Flavour fl) {
  return Kabbala(string("M_{"+fl.TexName()+"}"),fl.Yuk());
}

Kabbala Interaction_Model_sLepton_EW::K_yuk_sign(Flavour fl) {
  char hi[3];
  sprintf(hi,"%i",fl.MassSign());
  return Kabbala(string(hi),fl.MassSign());
}

Kabbala Interaction_Model_sLepton_EW::K_Z_H(short int i,short int j) {   
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("Z^{")+string(hi)+string(hj)+string("}_H"),
		 ComplexMatrixElement(string("Z_H"),i,j));
}     

Kabbala Interaction_Model_sLepton_EW::K_Z_R(short int i,short int j) {   
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("Z^{")+string(hi)+string(hj)+string("}_R"),
		 ComplexMatrixElement(string("Z_R"),i,j));
}  

Kabbala Interaction_Model_sLepton_EW::K_B_R(short int i) {   
  char hi[2];
  sprintf(hi,"%i",i);
  return Kabbala(string("B^{")+string(hi)+string("}_R"),
		 v1.Value() * ComplexMatrixElement(string("Z_R"),0,i) -
		 v2.Value() * ComplexMatrixElement(string("Z_R"),1,i) );
}  

Kabbala Interaction_Model_sLepton_EW::K_Z_PL(short int i,short int j)       
{   
  char hi[2];
  char hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("Z^\\p_{")+string(hi)+string(hj)+string("}"),
		 ComplexMatrixElement(string("Z^+"),i,j));
}  

Kabbala Interaction_Model_sLepton_EW::K_Z_MI(short int i,short int j)       
{   
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("Z^\\m_{")+string(hi)+string(hj)+string("}"),
		 ComplexMatrixElement(string("Z^-"),i,j));
}  

Kabbala Interaction_Model_sLepton_EW::K_Z_N(short int i,short int j)       
{   
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("Z^{")+string(hi)+string(hj)+string("}_N"),
		 ComplexMatrixElement(string("Z^N"),i,j));
}  
//we use transposed convention !!! 

Kabbala Interaction_Model_sLepton_EW::K_Z_N_com(short int i,short int j)       
{   
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);

  Complex exp_i = Complex(1.,0.);
  
  return Kabbala(string("Z^{")+string(hi)+string(hj)+string("}_N"),
		 exp_i*ComplexMatrixElement(string("Z^N"),i,j));
}  

Kabbala Interaction_Model_sLepton_EW::K_Z_N_com_conj(short int i,short int j)       
{   
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  Complex exp_i = Complex(1.,0.);
  
  return Kabbala(string("Z^{")+string(hi)+string(hj)+string("}_N"),
		 exp_i*ComplexMatrixElement(string("Z^N"),i,j));
}  


int Interaction_Model_sLepton_EW::gen_sLep(Flavour fl)
{
  int gen_sL;

  if (fl.Kfcode() == 71 || fl.Kfcode() == 74)
    gen_sL = 0;
  if (fl.Kfcode() == 72 || fl.Kfcode() == 75)
    gen_sL = 1;
  if (fl.Kfcode() == 73 || fl.Kfcode() == 76)
    gen_sL = 2;

  if (fl.Kfcode() == 81)
    gen_sL = 0;
  if (fl.Kfcode() == 82)
    gen_sL = 1;
  if (fl.Kfcode() == 83)
    gen_sL = 2;
  
  return gen_sL;
}

