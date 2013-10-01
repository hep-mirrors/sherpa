#include "MODEL/Interaction_Models/Interaction_Model_sLepton_EW.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include <stdio.h>


using namespace MODEL;
using namespace ATOOLS;
using namespace std;

Interaction_Model_sLepton_EW::Interaction_Model_sLepton_EW(MODEL::Model_Base * _model,
							   std::string _cplscheme,std::string _yukscheme) :
  Interaction_Model_Base("",_model,_cplscheme,_yukscheme)
{ 
  g1       = Kabbala(string("g_1"),
		     sqrt(4.*M_PI*ScalarFunction(string("alpha_QED"),rpa->gen.CplScale())));
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

  PL       = Kabbala(string("P_L"),1.);
  PR       = Kabbala(string("P_R"),1.);
  M_I      = Kabbala(string("i"),Complex(0.,1.));
  
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

void Interaction_Model_sLepton_EW::c_SSS(std::vector<Single_Vertex>& vertex,int& vanz)
{
  Kabbala kcpl0,kcpl1;
  int s_lep;
  
  //sneutrino - Higgs - sneutrino
  for (short int i=0;i<3;i++) {
    Flavour flav = Flavour((kf_code)(1000012+2*i));
    for (short int k=0;k<2;k++) {
      Flavour flh = Flavour((kf_code)(25+k*10));
      if (flh.IsOn() && flav.IsOn()) {
	vertex[vanz].in[0] = flav;
	vertex[vanz].in[1] = flh;
	vertex[vanz].in[2] = flav;

	kcpl0 = -M_I*g2*g2/(costW*costW*num_4)*K_B_R(k);
	kcpl1 = kcpl0;
	
	vertex[vanz].cpl[0]  = kcpl0; 
	vertex[vanz].cpl[1]  = kcpl1;
	vertex[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();

	
	vertex[vanz].Color.push_back(Color_Function(cf::None)); 
	
	
	vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("SSS",LF_Key()));
	
	vertex[vanz].on      = 1;
	vertex.push_back(Single_Vertex());vanz++;
      }
    }
  }

  //sneutrino - Hmin - slepton
 
  Flavour flHm = Flavour(kf_Hplus).Bar();
  if (flHm.IsOn()) {
    for (short int i=0;i<3;i++) {
      Flavour flav1 = Flavour((kf_code)(1000012+2*i));
      for (short int j=1;j<7;j++) {
	if (j<4) s_lep = 1000010 + 2*j - 1;
	else     s_lep = 2000010 + 2*j - 7;
	Flavour flav2 =Flavour((kf_code)(s_lep));
	if(flav1.IsOn() && flav2.IsOn()){
	  vertex[vanz].in[0] = flav1.Bar();
	  vertex[vanz].in[1] = flHm;
	  vertex[vanz].in[2] = flav2.Bar();
     
          Kabbala K_lI = Kabbala(string("\\frac{\\m M_{")+flav2.TexName()+string("}}{ v_1}\\sqrt{2}"),
			-K_yuk(Flavour((kf_code)(2*gen_sLep(flav2)+11))).Value()/v1.Value()*sqrt(2.));

	  kcpl0 = M_I*K_Z_Nu(gen_sLep(flav2),i)*
	    (-root2*g2*g2/num_4*(v1*K_Z_H(0,0)+v2*K_Z_H(1,0))*
	     K_Z_L(gen_sLep(flav2),j-1)-
	     K_Z_H(0,0)*(K_yuk(Flavour((kf_code)(2*gen_sLep(flav2)+11)))*K_lI*K_Z_L(gen_sLep(flav2),j-1)-
			 (K_l_S(gen_sLep(flav2),0)*K_Z_L(3,j-1)+
			  K_l_S(gen_sLep(flav2),1)*K_Z_L(4,j-1)+
			  K_l_S(gen_sLep(flav2),2)*K_Z_L(5,j-1)))
	     +K_Z_H(1,0)*(K_k_S(gen_sLep(flav2),0)*K_Z_L(3,j-1)+
			  K_k_S(gen_sLep(flav2),1)*K_Z_L(4,j-1)+
			  K_k_S(gen_sLep(flav2),2)*K_Z_L(5,j-1)-
			  K_lI*mu*K_Z_L(gen_sLep(flav2)+3,j-1)));
	 
	  kcpl1 = kcpl0;
	  
	  vertex[vanz].cpl[0]  = kcpl0; 
	  vertex[vanz].cpl[1]  = kcpl1;
	  vertex[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();
	
	  
	  vertex[vanz].Color.push_back(Color_Function(cf::None)); 
	      
	  
	  vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("SSS",LF_Key()));
	  
	  vertex[vanz].on      = 1;
	  vertex.push_back(Single_Vertex());vanz++;
	  
	}
      }
    }
  }
  
  //slepton - A0 - slepton
  Flavour flA0 = Flavour(kf_A0);
  if (flA0.IsOn()) {
    for (short int i=1;i<7;i++) {
      if (i<4) s_lep = 1000010 + 2*i - 1;
      else     s_lep = 2000010 + 2*i - 7;
      Flavour flav1 = Flavour((kf_code)(s_lep));
      for (short int j=i;j<7;j++) {
	if (j<4) s_lep = 1000010 + 2*j - 1;
	else     s_lep = 2000010 + 2*j - 7;
	Flavour flav2 =Flavour((kf_code)(s_lep));
	if(flav1.IsOn() && flav2.IsOn()){
	  vertex[vanz].in[0] = flav1;
	  vertex[vanz].in[1] = flA0;
	  vertex[vanz].in[2] = flav2;

	  Kabbala K_lI = Kabbala(string("\\frac{(\\m M_{")+flav1.TexName()+string("})}{ v_1}\\sqrt{2}"),
				 -K_yuk(Flavour((kf_code)(2*gen_sLep(flav1)+11))).Value()/(v1).Value()*sqrt(2.));
     
	  kcpl0 = -(K_l_S(gen_sLep(flav1),gen_sLep(flav2))*
		   (K_Z_L(gen_sLep(flav1),j-1)*K_Z_L(gen_sLep(flav2)+3,i-1)-
		    K_Z_L(gen_sLep(flav1),i-1)*K_Z_L(gen_sLep(flav2)+3,j-1))*K_Z_H(0,0)
		   + K_k_S(gen_sLep(flav1),gen_sLep(flav2))*
		   (K_Z_L(gen_sLep(flav1),j-1)*K_Z_L(gen_sLep(flav2)+3,i-1)-
		    K_Z_L(gen_sLep(flav1),i-1)*K_Z_L(gen_sLep(flav2)+3,j-1))*K_Z_H(1,0)
		   + K_lI*mu*(K_Z_L(gen_sLep(flav1),i-1)*K_Z_L(gen_sLep(flav1)+3,j-1)-
			      K_Z_L(gen_sLep(flav1),j-1)*K_Z_L(gen_sLep(flav1)+3,i-1))*
		   K_Z_H(1,0))*invroot2;
	  	  	      
	  kcpl1 = kcpl0;
	  
	  vertex[vanz].cpl[0]  = kcpl0; 
	  vertex[vanz].cpl[1]  = kcpl1;
	  vertex[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();
	  
	  
	  vertex[vanz].Color.push_back(Color_Function(cf::None)); 
	  
	  
	  vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("SSS",LF_Key()));
	  	  
	  vertex[vanz].on      = 1;
	  vertex.push_back(Single_Vertex());vanz++;
  
	}
      }
    }
  }
			       
  //slepton - h0/H0 - slepton
 
  for (short int k=0;k<2;k++) {
    Flavour flH = Flavour((kf_code)(25+k*10)); 
    for (short int i=1;i<7;i++) {
      if (i<4) s_lep = 1000010 + 2*i - 1;
      else     s_lep = 2000010 + 2*i - 7;
      Flavour flav1 = Flavour((kf_code)(s_lep));
      for (short int j=i;j<7;j++) {
	if (j<4) s_lep = 1000010 + 2*j - 1;
	else     s_lep = 2000010 + 2*j - 7;
	Flavour flav2 =Flavour((kf_code)(s_lep));
	if(flH.IsOn() && flav1.IsOn() && flav2.IsOn() && gen_sLep(flav1)==gen_sLep(flav2)){
	  vertex[vanz].in[0] = flav1;
	  vertex[vanz].in[1] = flH;
	  vertex[vanz].in[2] = flav2;

	  Kabbala help = K_zero;

	  Kabbala K_lI = Kabbala(string("\\frac{\\m M_{")+flav1.TexName()+string("}}{ v_1}*\\sqrt{2}"),
				 -K_yuk(Flavour((kf_code)(2*gen_sLep(flav1)+11))).Value()/(v1).Value()*sqrt(2.));
	 
	  Kabbala fac = Kabbala(string("\\frac{1-4sin^2\\theta_W}{2sin^2\\theta_W}"),
				(1.-4.*(sintW).Value()*(sintW).Value())/
				(2.*(sintW).Value()*(sintW).Value()));

	  if (i==j) {help = Kabbala(string("1"),1.);}
 
	  kcpl0 = M_I*(g1*g1/(costW*costW*num_2)*K_B_R(k)*
		       (help+fac*K_Z_L(gen_sLep(flav1),i-1)*K_Z_L(gen_sLep(flav1),j-1))
		       - K_lI*K_lI*v1*K_Z_R(0,k)*
		       (K_Z_L(gen_sLep(flav1),i-1)*K_Z_L(gen_sLep(flav1),j-1)+
		 	K_Z_L(gen_sLep(flav1)+3,i-1)*K_Z_L(gen_sLep(flav1)+3,j-1))
		       - K_Z_R(0,k)/root2*K_l_S(gen_sLep(flav1),gen_sLep(flav2))*
		       (K_Z_L(gen_sLep(flav1),j-1)*
			K_Z_L(gen_sLep(flav2)+3,i-1)+
			K_Z_L(gen_sLep(flav1),i-1)*
			K_Z_L(gen_sLep(flav2)+3,j-1))
		       + K_Z_R(1,k)/root2*K_k_S(gen_sLep(flav1),gen_sLep(flav2))*
		       (K_Z_L(gen_sLep(flav1),j-1)*
			K_Z_L(gen_sLep(flav2)+3,i-1)+
			K_Z_L(gen_sLep(flav1),i-1)*
			K_Z_L(gen_sLep(flav2)+3,j-1))
		       - K_lI*K_Z_R(1,k)*mu/root2*(K_Z_L(gen_sLep(flav1),i-1)*
						      K_Z_L(gen_sLep(flav1)+3,j-1)+
						      K_Z_L(gen_sLep(flav1),j-1)*
						      K_Z_L(gen_sLep(flav1)+3,i-1)));
						    

	  kcpl1 = kcpl0;

	  vertex[vanz].cpl[0]  = kcpl0;
	  vertex[vanz].cpl[1]  = kcpl1;
	  vertex[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();
	  	  
	  
	  vertex[vanz].Color.push_back(Color_Function(cf::None)); 
	  
	  
	  vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("SSS",LF_Key()));
	  
	  vertex[vanz].on      = 1;
	  vertex.push_back(Single_Vertex());vanz++;
  	}
      }
    }
  }
}

void Interaction_Model_sLepton_EW::c_SSV(std::vector<Single_Vertex>& vertex,int& vanz)
{
  Kabbala kcpl0,kcpl1;
  int s_lep;
  
  //slepton - Photon - slepton
  Flavour flPh = Flavour(kf_photon);
  if (flPh.IsOn()) {
    for (short int i=1;i<7;i++) {
      if (i<4) s_lep = 1000010 + 2*i - 1;
      else     s_lep = 2000010 + 2*i - 7;
      Flavour flav = Flavour((kf_code)(s_lep));
      Kabbala charge1 = Kabbala(string("Q_{")+flav.TexName()+string("}"),flav.Charge());
      if (flav.IsOn()) {
	vertex[vanz].in[0] = flav;
	vertex[vanz].in[1] = flPh;
	vertex[vanz].in[2] = flav;
	
	kcpl0 = -M_I*g1*charge1;
	kcpl1 = kcpl0;
	
	vertex[vanz].cpl[0]  = kcpl0; 
	vertex[vanz].cpl[1]  = kcpl1;
	vertex[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();

	
	vertex[vanz].Color.push_back(Color_Function(cf::None)); 
	
	
	vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("SSV",LF_Key()));
	vertex[vanz].Lorentz.back()->SetParticleArg(0,2,1);     	

	vertex[vanz].on      = 1;
	vertex.push_back(Single_Vertex());vanz++;
      }
    }
  }
 
 //slepton - Z - slepton
  Flavour flZ = Flavour(kf_Z);
  if (flZ.IsOn()) {
    for (short int i=1;i<7;i++) {
      if (i<4) s_lep = 1000010 + 2*i - 1;
      else     s_lep = 2000010 + 2*i - 7;
      Flavour flav1 = Flavour((kf_code)(s_lep));
      for (short int j=i;j<7;j++) {
	 if (j<4) s_lep = 1000010 + 2*j - 1;
	 else     s_lep = 2000010 + 2*j - 7;
	Flavour flav2 = Flavour((kf_code)(s_lep));
	if (flav1.IsOn() && flav2.IsOn() && gen_sLep(flav1)==gen_sLep(flav2)) {
	  
	  Kabbala help = K_zero;
	  if(i==j) help = sintW*sintW;
	  
	  vertex[vanz].in[0] = flav1;
	  vertex[vanz].in[1] = flZ;
	  vertex[vanz].in[2] = flav2;

	  kcpl0 = M_I*g2/costW*
	    ((K_Z_L(gen_sLep(flav2),j-1)*K_Z_L(gen_sLep(flav2),i-1))/num_2 - help);
	  kcpl1 = kcpl0;	  	  

	  vertex[vanz].cpl[0]  = kcpl0; 
	  vertex[vanz].cpl[1]  = kcpl1;
	  vertex[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();

	  
	  vertex[vanz].Color.push_back(Color_Function(cf::None)); 
	
	  
	  vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("SSV",LF_Key()));
	  vertex[vanz].Lorentz.back()->SetParticleArg(0,2,1);     	
	  
	  vertex[vanz].on      = 1;
	  vertex.push_back(Single_Vertex());vanz++;
	}
      }

    }

    //sneutrino - Z - sneutrino
    for (short int i=0;i<3;i++) {
      Flavour flav = Flavour((kf_code)(1000012+2*i));
      if (flav.IsOn()) {
 	  vertex[vanz].in[0] = flav;
	  vertex[vanz].in[1] = flZ;
	  vertex[vanz].in[2] = flav;

	  kcpl0 = -M_I*g2/(costW*num_2);
	  kcpl1 = kcpl0;
	  
	  vertex[vanz].cpl[0]  = kcpl0; 
	  vertex[vanz].cpl[1]  = kcpl1;
	  vertex[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();
	  
	  
	  vertex[vanz].Color.push_back(Color_Function(cf::None)); 
	  
	  
	  vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("SSV",LF_Key()));
	  vertex[vanz].Lorentz.back()->SetParticleArg(0,2,1);     	
	  
	  vertex[vanz].on      = 1;
	  vertex.push_back(Single_Vertex());vanz++;
      }
    }
    
  }
  
  //sneutrino - W - slepton
  Flavour flWplus = Flavour(kf_Wplus);
  if (flWplus.IsOn()) {
      for (short int i=0;i<3;i++) {
      Flavour flav1 = Flavour((kf_code)(1000012+2*i));
      for (short int j=1;j<7;j++) {
	if (j<4) s_lep = 1000010 + 2*j - 1;
	else     s_lep = 2000010 + 2*j - 7;
	Flavour flav2 =Flavour((kf_code)(s_lep));
	if(flav1.IsOn() && flav2.IsOn() && gen_sLep(flav1)==gen_sLep(flav2)) {
	  vertex[vanz].in[0] = flav1;
	  vertex[vanz].in[1] = flWplus;
	  vertex[vanz].in[2] = flav2;
	
	  kcpl0 = -M_I*g2/root2
	    *K_Z_Nu(gen_sLep(flav1),gen_sLep(flav1))
	    *K_Z_L(gen_sLep(flav1),j-1);
	  kcpl1 = kcpl0;
	  
	  vertex[vanz].cpl[0]  = kcpl0; 
	  vertex[vanz].cpl[1]  = kcpl1;
	  vertex[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();
	  
	  
	  vertex[vanz].Color.push_back(Color_Function(cf::None)); 
	  
	  
	  vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("SSV",LF_Key()));
	  vertex[vanz].Lorentz.back()->SetParticleArg(0,2,1);     	
	  
	  vertex[vanz].on      = 1;
	  vertex.push_back(Single_Vertex());vanz++;
  	}
      }
    }
  }
}

void Interaction_Model_sLepton_EW::c_SSSS(std::vector<Single_Vertex>& vertex,int& vanz)
{
  Flavour flHmin = Flavour(kf_Hplus).Bar();    
  Flavour flA0(kf_A0);    
  int s_lep;
  
  Kabbala kcpl0,kcpl1,help,K_lI, num_1;
  num_1    = Kabbala(string("1"),1.);    	

  for (short int i=1;i<7;++i) {
    if (i<4) s_lep = 1000010 + 2*i - 1;
    else     s_lep = 2000010 + 2*i - 7;
    Flavour flav1 = Flavour((kf_code)(s_lep));
      for (short int j=i;j<7;++j) {
	if (j<4) s_lep = 1000010 + 2*j - 1;
	else     s_lep = 2000010 + 2*j - 7;
	Flavour flav2 = Flavour((kf_code)(s_lep));
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
			   -K_yuk(Flavour((kf_code)(2*gen_sLep(flav1)+11))).Value()/(v1).Value()*sqrt(2.));
	    
	    
	    kcpl0 = M_I*(g1*g1/(costW*costW*num_2)*K_A_H(0,0)*
			 (help-(num_1+num_2*sintW*sintW)/(num_2*sintW*sintW)*
			  K_Z_L(gen_sLep(flav1),i-1)*K_Z_L(gen_sLep(flav1),j-1))-
			 K_lI*K_lI*K_Z_H(0,0)*K_Z_H(0,0)*
			 K_Z_L(gen_sLep(flav1)+3,i-1)*K_Z_L(gen_sLep(flav1)+3,j-1));
	    
	    kcpl1 = kcpl0;
	    
	    vertex[vanz].cpl[0]  = kcpl0; 
	    vertex[vanz].cpl[1]  = kcpl1;
	    vertex[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();
	    
	    
	    vertex[vanz].Color.push_back(Color_Function(cf::None)); 
	    
	    
	    vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("SSSS",LF_Key()));     
	    
	    vertex[vanz].on      = 1;
	    vertex[vanz].oew     = 2;
	    if (kcpl0.Value()!=Complex(0.,0.) && kcpl1.Value()!=Complex(0.,0.)) {vertex.push_back(Single_Vertex());vanz++;}
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
			   -K_yuk(Flavour((kf_code)(2*gen_sLep(flav1)+11))).Value()/(v1).Value()*sqrt(2.));
	    
	    
	    kcpl0 = M_I*(g1*g1/(costW*costW*num_2)*K_A_H(0,0)*
			 (help+(num_1-num_4*sintW*sintW)/(num_2*sintW*sintW)*
			  K_Z_L(gen_sLep(flav1),i-1)*K_Z_L(gen_sLep(flav1),j-1))-
			 K_lI*K_lI*K_Z_H(0,0)*K_Z_H(0,0)*
			 (K_Z_L(gen_sLep(flav1),i-1)*K_Z_L(gen_sLep(flav1),j-1) +
			 K_Z_L(gen_sLep(flav1)+3,i-1)*K_Z_L(gen_sLep(flav1)+3,j-1)));
	    
	    kcpl1 = kcpl0;
	    
	    vertex[vanz].cpl[0]  = kcpl0; 
	    vertex[vanz].cpl[1]  = kcpl1;
	    vertex[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();
	    
	    
	    vertex[vanz].Color.push_back(Color_Function(cf::None)); 
	    
	    
	    vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("SSSS",LF_Key()));     
	    
	    vertex[vanz].on      = 1;
	    if (kcpl0.Value()!=Complex(0.,0.) && kcpl1.Value()!=Complex(0.,0.)) {vertex.push_back(Single_Vertex());vanz++;}
	  }
	  //h0/H0 -> slepton - slepton - h0/H0 
	  for (int k=0;k<2;k++) {
	    Flavour flh1 = Flavour((kf_code)(25+k*10));
	    for (int l=k;l<2;l++) {
	      Flavour flh2 = Flavour((kf_code)(25+l*10));
	      if (flh1.IsOn() && flh2.IsOn()) {
		
		vertex[vanz].in[0] = flh1;
		vertex[vanz].in[1] = flav1.Bar();
		vertex[vanz].in[2] = flav2;
		vertex[vanz].in[3] = flh2;
		
		vertex[vanz].nleg  = 4;  
		
		help = K_zero;
		if(i==j) help = num_1;
		
		K_lI = Kabbala(string("\\frac{\\m M_{")+flav1.TexName()+string("}}{ v_1}*\\sqrt{2}"),
			       -K_yuk(Flavour((kf_code)(2*gen_sLep(flav1)+11))).Value()/(v1).Value()*sqrt(2.));
		
		
		kcpl0 = M_I*(g1*g1/(costW*costW*num_2)*K_A_R(k,l)*
			     (help+(num_1-num_4*sintW*sintW)/(num_2*sintW*sintW)*
			      K_Z_L(gen_sLep(flav1),i-1)*K_Z_L(gen_sLep(flav1),j-1))-
			     K_lI*K_lI*K_Z_R(0,k)*K_Z_R(0,l)*
			     (K_Z_L(gen_sLep(flav1),i-1)*K_Z_L(gen_sLep(flav1),j-1) +
			      K_Z_L(gen_sLep(flav1)+3,i-1)*K_Z_L(gen_sLep(flav1)+3,j-1)));
		
		kcpl1 = kcpl0;
		
		vertex[vanz].cpl[0]  = kcpl0; 
		vertex[vanz].cpl[1]  = kcpl1;
		vertex[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();
		
		
		vertex[vanz].Color.push_back(Color_Function(cf::None)); 
		
		
		vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("SSSS",LF_Key()));     
		
		vertex[vanz].on      = 1;
		vertex[vanz].oew     = 2;
		if (kcpl0.Value()!=Complex(0.,0.) && kcpl1.Value()!=Complex(0.,0.)) {vertex.push_back(Single_Vertex());vanz++;}
	      }
	    }
	  }
	}
      }
  }

  for (short int i=0;i<3;++i) {
    Flavour flsnu = Flavour((kf_code)(1000012+2*i));
    if (flsnu.IsOn()) {
      //h0/H0 -> snu - snub - h0/H0
        for (int k=0;k<2;k++) {
	    Flavour flh1 = Flavour((kf_code)(25+k*10));
	    for (int l=k;l<2;l++) {
	      Flavour flh2 = Flavour((kf_code)(25+l*10));
	      if (flh1.IsOn() && flh2.IsOn()) {
		
		vertex[vanz].in[0] = flh1;
		vertex[vanz].in[1] = flsnu.Bar();
		vertex[vanz].in[2] = flsnu;
		vertex[vanz].in[3] = flh2;
		
		vertex[vanz].nleg  = 4;  
		
		kcpl0 = -M_I*g2*g2/(costW*costW*num_4)*K_A_R(k,l);
		
		kcpl1 = kcpl0;
		
		vertex[vanz].cpl[0]  = kcpl0; 
		vertex[vanz].cpl[1]  = kcpl1;
		vertex[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();
		
		
		vertex[vanz].Color.push_back(Color_Function(cf::None)); 
		
		
		vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("SSSS",LF_Key()));     
		
		vertex[vanz].on      = 1;
		vertex[vanz].oew     = 2;
		vertex.push_back(Single_Vertex());vanz++;
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
	  
	  vertex[vanz].cpl[0]  = kcpl0; 
	  vertex[vanz].cpl[1]  = kcpl1;
	  vertex[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();
	  
	  
	  vertex[vanz].Color.push_back(Color_Function(cf::None)); 
	  
	  
	  vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("SSSS",LF_Key()));     
	  
	  vertex[vanz].on      = 1;
	  vertex[vanz].oew     = 2;
	  vertex.push_back(Single_Vertex());vanz++;
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
	  
	  
	  Flavour lepton = Flavour((kf_code)(11+i*2));
	  
	  Kabbala K_lI = Kabbala(string("\\frac{\\m M_{")+lepton.TexName()+string("}}{ v_1}\\sqrt{2}"),
				 -K_yuk(lepton).Value()/v1.Value()*sqrt(2.));
	  
	  
	  vertex[vanz].nleg  = 4;  
	  
	  kcpl0 = M_I*K_Z_Nu(i,i)*K_Z_Nu(i,i)*
	    (g1*g1*cot2TW/(num_2*sintW*costW)*K_A_H(0,0) - K_lI*K_lI*K_Z_H(0,0)*K_Z_H(0,0));
	  
	  kcpl1 = kcpl0;
	  
	  vertex[vanz].cpl[0]  = kcpl0; 
	  vertex[vanz].cpl[1]  = kcpl1;
	  vertex[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();
	  
	  
	  vertex[vanz].Color.push_back(Color_Function(cf::None)); 
	  
	  
	  vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("SSSS",LF_Key()));     
	  
	  vertex[vanz].on      = 1;
	  vertex[vanz].oew     = 2;
	  vertex.push_back(Single_Vertex());vanz++;
	} 
	if (flHmin.IsOn()) {
	  for (int j=1;j<7;j++) {
	    if (j<4) s_lep = 1000010 + 2*j - 1;
	    else     s_lep = 2000010 + 2*j - 7;
	    Flavour flslep = Flavour((kf_code)(s_lep));
	    if (gen_sLep(flsnu)==gen_sLep(flslep)) {
	      // H- -> sLepton - snub - h0/H0
	      for (int k=0;k<2;k++) {
		Flavour flh = Flavour((kf_code)(25+k*10));
		if (flh.IsOn()) {
		  vertex[vanz].in[0] = flHmin;
		  vertex[vanz].in[1] = flsnu.Bar();
		  vertex[vanz].in[2] = flslep;
		  vertex[vanz].in[3] = flh;
		  
		  Flavour lepton = Flavour((kf_code)(11+gen_sLep(flslep)*2));
		  
		  Kabbala K_lI = Kabbala(string("\\frac{\\m M_{")+lepton.TexName()+string("}}{ v_1}\\sqrt{2}"),
					 -K_yuk(lepton).Value()/v1.Value()*sqrt(2.));
		  
		  vertex[vanz].nleg  = 4;  
		  
		  kcpl0 = M_I/root2*K_Z_L(gen_sLep(flslep),j-1)*K_Z_Nu(gen_sLep(flslep),gen_sLep(flsnu))*
		    (-g2*g2/num_2*(K_Z_H(0,0)*K_Z_R(0,k) + K_Z_H(1,0)*K_Z_R(1,k)) 
		     + K_lI*K_lI*K_Z_H(0,0)*K_Z_R(0,k));
		  
		  kcpl1 = kcpl0;
	      
		  vertex[vanz].cpl[0]  = kcpl0; 
		  vertex[vanz].cpl[1]  = kcpl1;
		  vertex[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();
		  
		  
		  vertex[vanz].Color.push_back(Color_Function(cf::None)); 
		  
		  
		  vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("SSSS",LF_Key()));     
		  
		  vertex[vanz].on      = 1;
		  vertex[vanz].oew     = 2;
		  if (kcpl0.Value()!=Complex(0.,0.) && kcpl1.Value()!=Complex(0.,0.)) {vertex.push_back(Single_Vertex());vanz++;}
		}
	      }
	      // H- -> sLepton - snub - A0
	      if (flA0.IsOn()) {
		vertex[vanz].in[0] = flHmin;
		vertex[vanz].in[1] = flsnu.Bar();
		vertex[vanz].in[2] = flslep;
		vertex[vanz].in[3] = flA0;
		
		Flavour lepton = Flavour((kf_code)(11+gen_sLep(flslep)*2));
		
		Kabbala K_lI = Kabbala(string("\\frac{\\m M_{")+lepton.TexName()+string("}}{ v_1}\\sqrt{2}"),
				       -K_yuk(lepton).Value()/v1.Value()*sqrt(2.));
		
		vertex[vanz].nleg  = 4;  
		
		kcpl0 = -num_1/root2*K_Z_L(gen_sLep(flslep),j-1)*K_Z_Nu(gen_sLep(flslep),gen_sLep(flsnu))*
		  (g2*g2/num_2*K_A_H(0,0) - K_lI*K_lI*K_Z_H(0,0)*K_Z_H(0,0));
		
		kcpl1 = kcpl0;
		
		vertex[vanz].cpl[0]  = kcpl0; 
		vertex[vanz].cpl[1]  = kcpl1;
		vertex[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();
		
		
		vertex[vanz].Color.push_back(Color_Function(cf::None)); 
		
		
		vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("SSSS",LF_Key()));     
		
		vertex[vanz].on      = 1;
		vertex[vanz].oew     = 2;
		if (kcpl0.Value()!=Complex(0.,0.) && kcpl1.Value()!=Complex(0.,0.)) {vertex.push_back(Single_Vertex());vanz++;}
	      }
	    }
	  }
	}
    }
  }
  //snu - snu - snu - snu
  for (short int i=0;i<3;++i) {
    Flavour snu1 = Flavour((kf_code)(1000012+2*i));
    if (snu1.IsOn()) {
      for (short int j=i;j<3;++j) {
	Flavour snu2 = Flavour((kf_code)(1000012+2*j));
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
      	  
	  vertex[vanz].cpl[0]  = kcpl0; 
	  vertex[vanz].cpl[1]  = kcpl1;
	  vertex[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();
	  
	  
	  vertex[vanz].Color.push_back(Color_Function(cf::None)); 
	  
	  
	  vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("SSSS",LF_Key()));     
	  
	  vertex[vanz].on      = 1;
	  vertex[vanz].oew     = 2;
	  if (kcpl0.Value()!=Complex(0.,0.) && kcpl1.Value()!=Complex(0.,0.)) {vertex.push_back(Single_Vertex());vanz++;}
	
	  //snu - sLep - snu - sLep

	  for (int k=1;k<7;k++) {
	    if (k<4) s_lep = 1000010 + 2*k - 1;
	    else     s_lep = 2000010 + 2*k - 7;
	    Flavour slep1 = Flavour((kf_code)(s_lep));
	    if (slep1.IsOn()) {
	      for (int l=k;l<7;l++) {
		if (l<4) s_lep = 1000010 + 2*l - 1;
		else     s_lep = 2000010 + 2*l - 7;
		Flavour slep2 = Flavour((kf_code)(s_lep));
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
		      addendum += (g2*g2/num_2*K_Z_L(m,l-1)*K_Z_L(n,k-1) + 
				   K_l(m)*K_l(n)*K_Z_L(m+3,l-1)*K_Z_L(n+3,k-1))
			*K_Z_Nu(n,gen_sLep(snu1))*K_Z_Nu(m,gen_sLep(snu2));
		    }
		  }
		  
		  kcpl0 = -M_I*(-g1*g1/(num_2*costW*costW)*
				(help1 + (num_1-num_4*sintW*sintW)/(num_2*sintW*sintW)*
				 (K_Z_L(0,l-1)*K_Z_L(0,k-1)+
				  K_Z_L(1,l-1)*K_Z_L(1,k-1)+
				  K_Z_L(2,l-1)*K_Z_L(2,k-1)))*help2 + addendum);
				
		  kcpl1 = kcpl0;
		  
		  vertex[vanz].cpl[0]  = kcpl0; 
		  vertex[vanz].cpl[1]  = kcpl1;
		  vertex[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();
		  
		  
		  vertex[vanz].Color.push_back(Color_Function(cf::None)); 
		  
		  
		  vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("SSSS",LF_Key()));     
		  
		  if (kcpl0.Value()!=Complex(0.,0.) && kcpl1.Value()!=Complex(0.,0.)) vertex[vanz].on      = 1;
		  vertex[vanz].oew     = 2;
		  vertex.push_back(Single_Vertex());vanz++;
		}
	      }
	    }
	  }
	}
      }
    }
  }
  //slep - slep - slep - slep 
  for (int i=1;i<7;i++) {
    if (i<4) s_lep = 1000010 + 2*i - 1;
    else     s_lep = 2000010 + 2*i - 7;
    Flavour slepi = Flavour((kf_code)(s_lep));
    if (slepi.IsOn()) {
      for (int j=1;j<7;j++) {
	if (j<4) s_lep = 1000010 + 2*j - 1;
	else     s_lep = 2000010 + 2*j - 7;
	Flavour slepj = Flavour((kf_code)(s_lep));
	if (slepj.IsOn()) {
	  for (int k=1;k<7;k++) {
	    if (k<4) s_lep = 1000010 + 2*k - 1;
	    else     s_lep = 2000010 + 2*k - 7;
	    Flavour slepk = Flavour((kf_code)(s_lep));
	    if (slepk.IsOn()) {
	      for (int l=1;l<7;l++) {
		if (l<4) s_lep = 1000010 + 2*l - 1;
		else     s_lep = 2000010 + 2*l - 7;
		Flavour slepl = Flavour((kf_code)(s_lep));
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
		    ZLsum_il += K_Z_L(t,i-1)*K_Z_L(t,l-1);
		    ZLsum_ik += K_Z_L(t,i-1)*K_Z_L(t,k-1);
		    ZLsum_jl += K_Z_L(t,j-1)*K_Z_L(t,l-1);
		    ZLsum_jk += K_Z_L(t,j-1)*K_Z_L(t,k-1);
		  }
		  
		  Kabbala addendum = K_zero;
		  
		  for (int m=0;m<3;m++) {
		    for (int n=0;n<3;n++) {
		      addendum += K_l(m)*K_l(n)*
			(K_Z_L(m+3,i-1)*K_Z_L(n,j-1) + K_Z_L(m+3,j-1)*K_Z_L(n,i-1))*
			(K_Z_L(m,l-1)*K_Z_L(n+3,k-1) + K_Z_L(m,k-1)*K_Z_L(n+3,l-1)); 
		    }
		  }
		  
		  kcpl0 = -M_I*(g1*g1/(costW*costW)*(d_il*d_jk+d_ik*d_jl)
				- num_3*g1*g1/(num_2*costW*costW)*
				(d_il*ZLsum_jk+d_ik*ZLsum_jl+d_jl*ZLsum_ik+d_jk*ZLsum_il) +
				g2*g2*(num_1+num_2*num_4*sintW*sintW)/(num_4*costW*costW)*
				(ZLsum_il*ZLsum_jl + ZLsum_ik*ZLsum_jl) + 
				addendum);
		  
		  kcpl1 = kcpl0;
		  
		  vertex[vanz].cpl[0]  = kcpl0; 
		  vertex[vanz].cpl[1]  = kcpl1;
		  vertex[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();
		  
		  
		  vertex[vanz].Color.push_back(Color_Function(cf::None)); 
		  
		  
		  vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("SSSS",LF_Key()));     
		  
		  vertex[vanz].on      = 1;
		  vertex[vanz].oew     = 2;
		  if (kcpl0.Value()!=Complex(0.,0.) && kcpl1.Value()!=Complex(0.,0.)) {vertex.push_back(Single_Vertex());vanz++;}
		}
	      }
	    }
	  }
	}
      }
    }
  }
}


void Interaction_Model_sLepton_EW::c_SSVV(std::vector<Single_Vertex>& vertex,int& vanz)
{
  Flavour flavWm = Flavour(kf_Wplus).Bar();
  Flavour flavZ(kf_Z);
  Flavour flavPhoton(kf_photon);
  Kabbala kcpl0,kcpl1,help,num_1;
  num_1    = Kabbala(string("1"),1.);    	
  int s_lep;

  for (short int i=1;i<7;i++) {
    if (i<4) s_lep = 1000010 + 2*i - 1;
    else     s_lep = 2000010 + 2*i - 7;
    Flavour flav = Flavour((kf_code)(s_lep));
    if(flav.IsOn() && flavPhoton.IsOn()){
      
      // P - L - L - P  
      vertex[vanz].in[0] = flavPhoton;
      vertex[vanz].in[1] = flav.Bar();
      vertex[vanz].in[2] = flav;
      vertex[vanz].in[3] = flavPhoton;
      
      vertex[vanz].nleg     = 4;
      
      kcpl0 = num_2*M_I*g1*g1;
      kcpl1 = kcpl0;
      
      vertex[vanz].cpl[0]  = kcpl0;
      vertex[vanz].cpl[1]  = kcpl1;
      vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	  
      
      vertex[vanz].Color.push_back(Color_Function(cf::None));     
      
      
      vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("VVSS",LF_Key()));     
      vertex[vanz].Lorentz.back()->SetParticleArg(0,3);     
      
      vertex[vanz].on      = 1;
      vertex[vanz].oew     = 2;
      vertex.push_back(Single_Vertex());vanz++;
    }
  }
  for (short int i=1;i<7;i++) {
    if (i<4) s_lep = 1000010 + 2*i - 1;
    else     s_lep = 2000010 + 2*i - 7;
    Flavour flav1 = Flavour((kf_code)(s_lep));
    for (short int j=1;j<7;j++) {
      if (j<4) s_lep = 1000010 + 2*j - 1;
      else     s_lep = 2000010 + 2*j - 7;
      Flavour flav2 = Flavour((kf_code)(s_lep));
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
	    
	    kcpl0 = M_I*g1*g2/costW*((K_Z_L(0,i-1)*K_Z_L(0,j-1) + 
				      K_Z_L(1,i-1)*K_Z_L(1,j-1) +
				      K_Z_L(2,i-1)*K_Z_L(2,j-1))
				      - num_2*help*sintW*sintW);
	    
	    kcpl1 = kcpl0;
	    
	    vertex[vanz].cpl[0]  = kcpl0;
	    vertex[vanz].cpl[1]  = kcpl1;
	    vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	    
	    
	    vertex[vanz].Color.push_back(Color_Function(cf::None));     
	    
	    
	    vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("VVSS",LF_Key()));     
	    vertex[vanz].Lorentz.back()->SetParticleArg(0,3);     
	    
	    vertex[vanz].on      = 1;
	    vertex[vanz].oew     = 2;
	    vertex.push_back(Single_Vertex());vanz++;
	  }			   
	  // Z - L - L - Z
	    vertex[vanz].in[0] = flavZ;
	    vertex[vanz].in[1] = flav1.Bar();
	    vertex[vanz].in[2] = flav2;
	    vertex[vanz].in[3] = flavZ;
	    
	    vertex[vanz].nleg     = 4;
	    
	    kcpl0 = num_2*M_I*g1*g1/(costW*costW)*(help*sintW*sintW +
						   (num_1 - num_4*sintW*sintW)/(num_4*sintW*sintW)*
						   (K_Z_L(0,i-1)*K_Z_L(0,j-1) +
						    K_Z_L(1,i-1)*K_Z_L(1,j-1) +
						    K_Z_L(2,i-1)*K_Z_L(2,j-1)));
	    
	    kcpl1 = kcpl0;
	    
	    vertex[vanz].cpl[0]  = kcpl0;
	    vertex[vanz].cpl[1]  = kcpl1;
	    vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	    
	    
	    vertex[vanz].Color.push_back(Color_Function(cf::None));     
	    
	    
	    vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("VVSS",LF_Key()));     
	    vertex[vanz].Lorentz.back()->SetParticleArg(0,3);     
	    
	    vertex[vanz].on      = 1;
	    vertex[vanz].oew     = 2;
	    vertex.push_back(Single_Vertex());vanz++;
	}
      }
    }
  }
  
  for (short int i=0;i<3;i++) {
    Flavour flav = Flavour((kf_code)(1000012+2*i));
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
	
	vertex[vanz].cpl[0]  = kcpl0;
	vertex[vanz].cpl[1]  = kcpl1;
	vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	
	
	vertex[vanz].Color.push_back(Color_Function(cf::None));     
	
	
	vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("VVSS",LF_Key()));
	vertex[vanz].Lorentz.back()->SetParticleArg(0,3);     
	
	vertex[vanz].on      = 1;
	vertex[vanz].oew     = 2;
	vertex.push_back(Single_Vertex());vanz++;
      }
      if (flavWm.IsOn()) {
	// W - snu - snu - W 
	vertex[vanz].in[0] = flavWm;
	vertex[vanz].in[1] = flav.Bar();
	vertex[vanz].in[2] = flav;
	vertex[vanz].in[3] = flavWm;
	
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
	vertex[vanz].oew     = 2;
	vertex.push_back(Single_Vertex());vanz++;
      }
    }
  }
  for (short int i=1;i<7;i++) {
    if (i<4) s_lep = 1000010 + 2*i - 1;
    else     s_lep = 2000010 + 2*i - 7;
    Flavour flav1 = Flavour((kf_code)(s_lep));
    for (short int j=i;j<7;j++) {
      if (j<4) s_lep = 1000010 + 2*j - 1;
      else     s_lep = 2000010 + 2*j - 7;
      Flavour flav2 =Flavour((kf_code)(s_lep));
      if(flav1.IsOn() && flav2.IsOn() && gen_sLep(flav1)==gen_sLep(flav2)){
	
	// W - L - L - W  
	if (flavWm.IsOn() && gen_sLep(flav1)==gen_sLep(flav2)) {
	  
	  vertex[vanz].in[0] = flavWm;
	  vertex[vanz].in[1] = flav1.Bar();
	  vertex[vanz].in[2] = flav2;
	  vertex[vanz].in[3] = flavWm;
	  
	  vertex[vanz].nleg     = 4;
	  
	  kcpl0 = M_I*g2*g2/num_2*K_Z_L(gen_sLep(flav1),i-1)*K_Z_L(gen_sLep(flav1),j-1);
	  kcpl1 = kcpl0;
	  
	  vertex[vanz].cpl[0]  = kcpl0;
	  vertex[vanz].cpl[1]  = kcpl1;
	  vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	  
	  
	  vertex[vanz].Color.push_back(Color_Function(cf::None));     
	  
	  
	  vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("VVSS",LF_Key()));
	  vertex[vanz].Lorentz.back()->SetParticleArg(0,3);     
	  
	  vertex[vanz].on      = 1;
	  vertex[vanz].oew     = 2;
	  vertex.push_back(Single_Vertex());vanz++;
	}
      }
    }
  }
  
  // W - snu - L - P/Z  
  for (short int i=0;i<3;i++) {
    Flavour flav1 = Flavour((kf_code)(1000012+2*i));
    for (short int j=1;j<7;j++) {
      if (j<4) s_lep = 1000010 + 2*j - 1;
      else     s_lep = 2000010 + 2*j - 7;
      Flavour flav2 =Flavour((kf_code)(s_lep));
      if(flav1.IsOn() && flav2.IsOn() && gen_sLep(flav1)==gen_sLep(flav2)){
	
	// W - snu - L - P  
	if (flavWm.IsOn()) {
	  if (flavPhoton.IsOn()) {
	    
	    vertex[vanz].in[0] = flavWm;
	    vertex[vanz].in[1] = flav1.Bar();
	    vertex[vanz].in[2] = flav2;
	    vertex[vanz].in[3] = flavPhoton;
	    
	    vertex[vanz].nleg     = 4;
	    
	    kcpl0 = -M_I*g1*g2/root2*(K_Z_Nu(gen_sLep(flav1),i)*
				      K_Z_L(gen_sLep(flav1),j-1));
	    kcpl1 = kcpl0;
	    
	    vertex[vanz].cpl[0]  = kcpl0;
	    vertex[vanz].cpl[1]  = kcpl1;
	    vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	    
	    
	    vertex[vanz].Color.push_back(Color_Function(cf::None));     
	    
	    
	    vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("VVSS",LF_Key()));
	    vertex[vanz].Lorentz.back()->SetParticleArg(0,3);     
	    
	    vertex[vanz].on      = 1;
	    vertex[vanz].oew     = 2;
	    vertex.push_back(Single_Vertex());vanz++;
	  }
	  
	  if (flavZ.IsOn()) {
	    // W - snu - L - Z  
	    
	    vertex[vanz].in[0] = flavWm;
	    vertex[vanz].in[1] = flav1.Bar();
	    vertex[vanz].in[2] = flav2;
	    vertex[vanz].in[3] = flavZ;
	    
	    vertex[vanz].nleg     = 4;
	    
	    kcpl0 = M_I*g1*g1/(root2*costW)*(K_Z_Nu(gen_sLep(flav1),i)*
					     K_Z_L(gen_sLep(flav1),j-1));
	    kcpl1 = kcpl0;
	    
	    vertex[vanz].cpl[0]  = kcpl0;
	    vertex[vanz].cpl[1]  = kcpl1;
	    vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	    
	    
	    vertex[vanz].Color.push_back(Color_Function(cf::None));     
	    
	    
	    vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("VVSS",LF_Key()));
	    vertex[vanz].Lorentz.back()->SetParticleArg(0,3);     
	    
	    vertex[vanz].on      = 1;
	    vertex[vanz].oew     = 2;
	    vertex.push_back(Single_Vertex());vanz++;
	  }
	}
      }
    }
  }
}

void Interaction_Model_sLepton_EW::c_FFS(std::vector<Single_Vertex>& vertex,int& vanz){

  Kabbala kcpl0,kcpl1;
  int n_ino,c_ino,s_lep;
  
  Kabbala K_l1 = Kabbala(string("l^1"),
			 -K_yuk(Flavour((kf_code)(11))).Value()/v1.Value()*sqrt(2.));

  Kabbala K_l2 = Kabbala(string("l^2"),
			 -K_yuk(Flavour((kf_code)(13))).Value()/v1.Value()*sqrt(2.));
  
  Kabbala K_l3 = Kabbala(string("l^3"),
			 -K_yuk(Flavour((kf_code)(15))).Value()/v1.Value()*sqrt(2.));
  
  //neutrino - sneutrino - neutralino
  
  for (short int i=12;i<17;i+=2) {
    Flavour flav1 = Flavour((kf_code)(i));
    for (short int j=0;j<4;j++) {
      if (j<2)  n_ino=1000022+j;
      if (j==2) n_ino=1000025;
      if (j==3) n_ino=1000035;
      Flavour flav2 = Flavour((kf_code)(n_ino));
      for (short int k=0;k<3;k++) {
	Flavour flav3 = Flavour((kf_code)(1000012+2*k));
	if (flav1.IsOn() && flav2.IsOn() && flav3.IsOn()) {
	  vertex[vanz].in[0] = flav1;
	  vertex[vanz].in[1] = flav3;
	  vertex[vanz].in[2] = flav2;
	  
	  kcpl0 = K_zero;
	  kcpl1 = M_I*g2/(costW*root2)*K_Z_Nu((i-12)/2,k)*
	    (K_Z_N(0,j)*sintW-K_Z_N(1,j)*costW);
	  
	  vertex[vanz].cpl[0] = kcpl0;
	  vertex[vanz].cpl[1] = kcpl1;
	  vertex[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();
	
	  
	  vertex[vanz].Color.push_back(Color_Function(cf::None)); 
	  
	  
	  vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("FFS",LF_Key()));
	  
	  vertex[vanz].on      = 1;
	  vertex.push_back(Single_Vertex());vanz++;
	}
      }
    }
  }  
    
  //neutrino - slepton - Chargino  

  for (short int i=12;i<17;i+=2) {
    Flavour flav1 = Flavour((kf_code)(i));
    
    Kabbala K_lI = Kabbala(string("\\frac{\\m M_{")+flav1.TexName()+string("}}{ v_1}\\sqrt{2}"),
			   -K_yuk(Flavour((kf_code)(i-1))).Value()/v1.Value()*sqrt(2.));
   
    for (short int j=0;j<2;j++) {
      if (j==0) c_ino=1000024;
      else      c_ino=1000037;
      Flavour flav2 = Flavour((kf_code)(c_ino));
      for (short int k=1;k<7;k++) {
	if (k<4) s_lep = 1000010 + 2*k - 1;
	else     s_lep = 2000010 + 2*k - 7;
	Flavour flav3 = Flavour((kf_code)(s_lep));
	if (flav1.IsOn() && flav2.IsOn() && flav3.IsOn()) {
	  vertex[vanz].in[0] = flav1;
	  vertex[vanz].in[1] = flav3;
	  vertex[vanz].in[2] = flav2;
	  
	  kcpl0 = K_zero;
	  kcpl1 = -M_I*(K_Z_MI(0,j)*K_Z_L((i-12)/2,k-1)*g2+
			K_Z_MI(1,j)*K_Z_L((i-12)/2+3,k-1)*K_lI);
	  
	  vertex[vanz].cpl[0] = kcpl0;
	  vertex[vanz].cpl[1] = kcpl1;
	  vertex[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();
	  
	  
	  vertex[vanz].Color.push_back(Color_Function(cf::None)); 
	  
	  
	  vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("FFS",LF_Key()));
	  
	  vertex[vanz].on     = 1;
	  vertex.push_back(Single_Vertex());vanz++;
	}
      }
    }
  }

  //lepton - Chargino - sneutrino  
 
  for (short int i=11;i<16;i+=2) {
    Flavour flav1 = Flavour((kf_code)(i));
  
    Kabbala K_lI = Kabbala(string("\\frac{(\\m M_{")+flav1.TexName()+string(")}}{ v_1}\\sqrt{2}"),
			   -K_yuk(flav1).Value()/v1.Value()*sqrt(2.));
       
    for (short int j=0;j<2;j++) {
      if (j==0) c_ino=1000024;
      else      c_ino=1000037;
      Flavour flav2 = Flavour((kf_code)(c_ino));
      for (short int k=0;k<3;k++) {
	Flavour flav3 = Flavour((kf_code)(1000012+2*k));
	if (flav1.IsOn() && flav2.IsOn() && flav3.IsOn()) {
	  vertex[vanz].in[0] = flav1;
	  vertex[vanz].in[1] = flav3;
	  vertex[vanz].in[2] = flav2.Bar();
	  
	  kcpl0 = -M_I*K_lI*K_Z_MI(1,j)*K_Z_Nu((i-11)/2,k);
	  kcpl1 = -M_I*g2*K_Z_PL(0,j)*K_Z_Nu((i-11)/2,k); 
	  
	  vertex[vanz].cpl[0] = kcpl0;
	  vertex[vanz].cpl[1] = kcpl1;
	  vertex[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();

	  
	  vertex[vanz].Color.push_back(Color_Function(cf::None)); 
	  
	  
	  vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("FFS",LF_Key()));
	  
	  vertex[vanz].on      = 1;
	  vertex.push_back(Single_Vertex());vanz++;
	}
      }
    }
  }

  //lepton - slepton - Neutralino  
 
  for (short int i=11;i<16;i+=2) {
    Flavour flav1 = Flavour((kf_code)(i));
    for (short int j=0;j<4;j++) {
      if (j<2)  n_ino=1000022+j;
      if (j==2) n_ino=1000025;
      if (j==3) n_ino=1000035;
      Flavour flav2 = Flavour((kf_code)(n_ino));
      for (short int k=1;k<7;k++) {
	if (k<4) s_lep = 1000010 + 2*k - 1;
	else     s_lep = 2000010 + 2*k - 7;
	Flavour flav3 = Flavour((kf_code)(s_lep));
	if (flav1.IsOn() && flav2.IsOn() && flav3.IsOn() && (i-11)/2==gen_sLep(flav3)) {
	  
	  //fixed + save
	  vertex[vanz].in[0] = flav1;
	  vertex[vanz].in[1] = flav3;
	  vertex[vanz].in[2] = flav2;

	  Kabbala K_lI = Kabbala(string("\\frac{\\m M_{")+flav1.TexName()+string("}}{ v_1}\\sqrt{2}"),
				     -K_yuk(flav1).Value()/v1.Value()*sqrt(2.));
	  
	  kcpl0 = M_I*(-g1/costW*root2*K_Z_L((i-11)/2+3,k-1)*K_Z_N_com_conj(0,j)
		       + K_lI*K_Z_L((i-11)/2,k-1)*K_Z_N_com_conj(2,j));
	  
	  kcpl1 = M_I*(g2/(costW*root2)*K_Z_L((i-11)/2,k-1)*
		       (K_Z_N_com(0,j)*sintW + K_Z_N_com(1,j)*costW) 
		       + K_lI*K_Z_L((i-11)/2+3,k-1)*K_Z_N_com(2,j));
	  
	  vertex[vanz].cpl[0] = kcpl0;
	  vertex[vanz].cpl[1] = kcpl1;
	  vertex[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();
	  
	  
	  vertex[vanz].Color.push_back(Color_Function(cf::None)); 
	  
	  
	  vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("FFS",LF_Key()));
	  
	  vertex[vanz].on     = 1;
	  vertex.push_back(Single_Vertex());vanz++;
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
		 -K_yuk(Flavour((kf_code)(2*i+11))).Value()/v1.Value()*sqrt(2.));

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
  if(ScalarNumber(std::string("WidthScheme"))==0){
    return Kabbala(string("M_{"+fl.TexName()+"}"),fl.Yuk());
  }else{
    return Kabbala(string("M_{"+fl.TexName()+"}"),
                   sqrt(sqr(fl.Yuk())-Complex(0.,1.)*fl.Width()*fl.Yuk()));
  }
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

  if (fl.Kfcode() == 1000011 || fl.Kfcode() == 2000011)
    gen_sL = 0;
  if (fl.Kfcode() == 1000013 || fl.Kfcode() == 2000013)
    gen_sL = 1;
  if (fl.Kfcode() == 1000015 || fl.Kfcode() == 2000015)
    gen_sL = 2;

  if (fl.Kfcode() == 1000012)
    gen_sL = 0;
  if (fl.Kfcode() == 1000014)
    gen_sL = 1;
  if (fl.Kfcode() == 1000016)
    gen_sL = 2;
  
  return gen_sL;
}

