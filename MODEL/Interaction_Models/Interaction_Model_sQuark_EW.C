#include "MODEL/Interaction_Models/Interaction_Model_sQuark_EW.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include <stdio.h>


using namespace MODEL;
using namespace ATOOLS;
using namespace std;

Interaction_Model_sQuark_EW::Interaction_Model_sQuark_EW(MODEL::Model_Base * _model,
							 std::string _cplscheme,std::string _yukscheme) :
  Interaction_Model_Base("",_model,_cplscheme,_yukscheme)
{ 
  double scale = rpa.gen.CplScale();

  g1       = Kabbala(string("g_1"),
		     sqrt(4.*M_PI*ScalarFunction(string("alpha_QED"),scale)));
  g2       = Kabbala(string("g_1/\\sin\\theta_W"), 
		     g1.Value()/sqrt(ScalarConstant(string("sin2_thetaW"))));
  g3    = Kabbala(string("g_3"),
		  sqrt(4.*M_PI*ScalarFunction(std::string("alpha_S"),scale)));
  sintW    = Kabbala(string("\\sin\\theta_W"),
		     sqrt(ScalarConstant(string("sin2_thetaW"))));
  costW    = Kabbala(string("\\cos\\theta_W"),
		     sqrt(1.-ScalarConstant(string("sin2_thetaW"))));
  PL       = Kabbala(string("P_L"),1.);
  PR       = Kabbala(string("P_R"),1.);
  M_I      = Kabbala(string("i"),Complex(0.,1.));
  vev      = Kabbala(string("v_{EW}"),ScalarConstant(string("vev")));
  v1       = Kabbala(string("v_1"),
		     ScalarConstant(string("vev")) *
		     sqrt(1./(1.+sqr(ScalarConstant(string("tan(beta)"))))));
  v2       = Kabbala(string("v_2"),
		     ScalarConstant(string("vev")) *
		     ScalarConstant(string("tan(beta)")) *
		     sqrt(1./(1.+sqr(ScalarConstant(string("tan(beta)"))))));
  mu       = Kabbala(string("h"),ScalarConstant(string("mu")));
  conj_mu  = Kabbala(string("h"),ScalarConstant(string("mu")));
  K_zero   = Kabbala(string("zero"),0.);
  num_1    = Kabbala(string("1"),1.);    	
  num_2    = Kabbala(string("2"),2.);    	
  num_3    = Kabbala(string("3"),3.);    	
  num_4    = Kabbala(string("4"),4.);    		
  num_6    = Kabbala(string("6"),6.);
  root2    = Kabbala(string("\\sqrt{2}"),sqrt(2.));
  invroot2 = Kabbala(string("1/\\sqrt{2}"),sqrt(.5));
}


void Interaction_Model_sQuark_EW::c_FFS(std::vector<Single_Vertex>& vertex,int& vanz) {  
  Kabbala kcpl0,kcpl1; 
  int n_ino,c_ino,s_qu;
  
  //quark - squark - neutralino
  for (short int j=0;j<4;j++) {
    if (j<2)  n_ino=1000022+j; 
    if (j==2) n_ino=1000025; 
    if (j==3) n_ino=1000035; 
    Flavour flneu = Flavour((kf_code)(n_ino));
    if (flneu.IsOn()) {
      //uptypes 
      for (short int k=2;k<7;k+=2) {
	Flavour flav1 = Flavour((kf_code)(k));
	for (short int i=1;i<7;i++) {
	  if (i<4) s_qu=1000000 + 2*i;
	  else     s_qu=2000000 + 2*i - 6;
	  Flavour flav2 = Flavour((kf_code)(s_qu));
	  if (flav1.IsOn() && flav2.IsOn() && (k/2-1)==gen_sUp(flav2)) {
	    vertex[vanz].in[0] = flav1;
	    vertex[vanz].in[1] = flav2;
	    vertex[vanz].in[2] = flneu;
	    
	    Kabbala K_uI = Kabbala(string("\\frac{\\m M_{")+flav1.TexName()+string("}}{ v_2}\\sqrt{2}"),
				     flav1.Yuk()/v2.Value()*sqrt(2.));
	  
	    kcpl0 = M_I*((g1*root2*num_2)/(costW*num_3)*K_Z_U((k-2)/2+3,i-1)*K_Z_N(0,j)+
			 -(K_uI*K_Z_U(gen_sUp(flav2),i-1))*K_Z_N(3,j)); 
	    
	    kcpl1 = M_I*(-g2/(costW*root2)*K_Z_U((k-2)/2,i-1)*
			 (K_Z_N(0,j)*(sintW/num_3) + K_Z_N(1,j)*costW)
 			 -(K_uI*K_Z_U(gen_sUp(flav2)+3,i-1))*K_Z_N(3,j));
	    
	    vertex[vanz].cpl[0] = kcpl0;
	    vertex[vanz].cpl[1] = kcpl1;
	    vertex[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();

	    
	    vertex[vanz].Color.push_back(Color_Function(cf::D));     
	    vertex[vanz].Color.back().SetParticleArg(0,1);     
	    vertex[vanz].Color.back().SetStringArg('0','1');     
	  
	    
	    vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("FFS",LF_Key()));
	    
	    vertex[vanz].on      = 1;
	    vertex.push_back(Single_Vertex());vanz++;
	  }
	}
      }
      //downtypes 
      for (short int k=1;k<6;k+=2) {
	Flavour flav1 = Flavour((kf_code)(k));
	for (short int i=1;i<7;i++) {
	  if (i<4) s_qu = 1000000 + 2*i -1;
	  else     s_qu = 2000000 + 2*i -7;
	  Flavour flav2 = Flavour((kf_code)(s_qu));
	  if (flav1.IsOn() && flav2.IsOn() && ((k-1)/2)==gen_sDown(flav2)) {
	    vertex[vanz].in[0] = flav1;
	    vertex[vanz].in[1] = flav2;
	    vertex[vanz].in[2] = flneu;
	    
	    Kabbala K_dI = Kabbala(string("d^I"),-flav1.Yuk()*sqrt(2.)/v1.Value());
	    
	    kcpl0 = M_I*((-(g1*root2)/
			  (costW*num_3))*K_Z_D((k-1)/2+3,i-1)*K_Z_N(0,j)+
			 (K_dI*K_Z_D(gen_sDown(flav2),i-1)*K_Z_N(2,j)));
	    
	    kcpl1 = M_I*(-g2/(costW*root2)*K_Z_D((k-1)/2,i-1)*
			 (K_Z_N(0,j)*(sintW/num_3)-K_Z_N(1,j)*costW)+
			 (K_dI*K_Z_D(gen_sDown(flav2)+3,i-1))*K_Z_N(2,j));
	    
	    vertex[vanz].cpl[0] = kcpl0;
	    vertex[vanz].cpl[1] = kcpl1;
	    vertex[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();
	    
	    
	    vertex[vanz].Color.push_back(Color_Function(cf::D));     
	    vertex[vanz].Color.back().SetParticleArg(0,1);     
	    vertex[vanz].Color.back().SetStringArg('0','1');     
	    
	    
	    vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("FFS",LF_Key()));
	    
	    vertex[vanz].on     = 1;
	    vertex.push_back(Single_Vertex());vanz++;
	  }
	}
      }
    }
  }
  
  //d-quark - Chargino - sup
  for (short int i=1;i<6;i+=2) {
    Flavour flav1 = Flavour((kf_code)(i));
    Kabbala K_dI = Kabbala(string("d^I"),
			   -flav1.Yuk()/v1.Value()*sqrt(2.));
    for (short int j=0;j<2;j++) {
      if (j==0) c_ino=1000024;
      else      c_ino=1000037;
      Flavour flav2 = Flavour((kf_code)(c_ino));
      for (short int k=1;k<7;k++) {
	if (k<4) s_qu=1000000 + 2*k;
	else     s_qu=2000000 + 2*k - 6;
	Flavour flav3 = Flavour((kf_code)(s_qu));
	if (flav1.IsOn() && flav2.IsOn() && flav3.IsOn()) {
	  vertex[vanz].in[0] = flav1;
	  vertex[vanz].in[1] = flav3;
	  vertex[vanz].in[2] = flav2.Bar();
	  
	  Kabbala K_uI = Kabbala(string("u^I"),Flavour((kf_code)(2*gen_sUp(flav3)+2)).Yuk()*
				 sqrt(2.)/v2.Value());

	  kcpl0 = -M_I*K_dI*K_Z_MI(1,j)*K_Z_U(gen_sUp(flav3),k-1)*
	    K_CKM(gen_sUp(flav3),(i-1)/2);
	  
	  kcpl1 = M_I*(-g2*K_Z_PL(0,j)*K_Z_U(gen_sUp(flav3),k-1)+
		       K_uI*K_Z_PL(1,j)*K_Z_U(gen_sUp(flav3)+3,k-1))*
		       K_CKM(gen_sUp(flav3),(i-1)/2);

	  vertex[vanz].cpl[0] = kcpl0;
	  vertex[vanz].cpl[1] = kcpl1;
	  vertex[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();
	  
	  
	  vertex[vanz].Color.push_back(Color_Function(cf::D));     
	  vertex[vanz].Color.back().SetParticleArg(0,1);     
	  vertex[vanz].Color.back().SetStringArg('0','1');     
	  
	  
	  vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("FFS",LF_Key()));
	  
	  vertex[vanz].on      = 1;
	  vertex.push_back(Single_Vertex());vanz++;
	}
     }
    }
  }
  //u-quark - Chargino - sdown
  for (short int i=2;i<7;i+=2) {
    Flavour flav1 = Flavour((kf_code)(i));
    Kabbala K_uJ = Kabbala(string("u^J"),flav1.Yuk()*sqrt(2.)/v2.Value());
    for (short int j=0;j<2;j++) {
      if (j==0) c_ino=1000024;
      else      c_ino=1000037;
      Flavour flav2 = Flavour((kf_code)(c_ino));
      for (short int k=1;k<7;k++) {
	if (k<4) s_qu=1000000 + 2*k - 1;
	else     s_qu=2000000 + 2*k - 7;
	Flavour flav3 = Flavour((kf_code)(s_qu));
	if (flav1.IsOn() && flav2.IsOn() && flav3.IsOn()) {
	  vertex[vanz].in[0] = flav1;
	  vertex[vanz].in[1] = flav3;
	  vertex[vanz].in[2] = flav2;
	  	  
	  Kabbala K_dI = Kabbala(string("d^I"),
				 -Flavour((kf_code)(2*gen_sDown(flav3)+1)).Yuk()*sqrt(2.)/
				 v1.Value());
	  
	  kcpl0 = M_I*K_uJ*K_Z_D(gen_sDown(flav3),k-1)*K_Z_PL(1,j)*
	    conj_K_CKM((i-2)/2,gen_sDown(flav3));
	  
	  kcpl1 = -M_I*(g2*K_Z_D(gen_sDown(flav3),k-1)*K_Z_MI(0,j)+
			K_dI*K_Z_D(gen_sDown(flav3)+3,k-1)*K_Z_MI(1,j))*
			conj_K_CKM((i-2)/2,gen_sDown(flav3));
			
	  vertex[vanz].cpl[0] = kcpl0;
	  vertex[vanz].cpl[1] = kcpl1;
	  vertex[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();
	  
	  
	  vertex[vanz].Color.push_back(Color_Function(cf::D));     
	  vertex[vanz].Color.back().SetParticleArg(0,1);     
	  vertex[vanz].Color.back().SetStringArg('0','1');     
	  
	  
	  vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("FFS",LF_Key()));
	  
	  vertex[vanz].on     = 1;
	  vertex.push_back(Single_Vertex());vanz++;
	}
      }
    } 
  }
}

void Interaction_Model_sQuark_EW::c_SSV(std::vector<Single_Vertex>& vertex,int& vanz)
{
  Kabbala kcpl0,kcpl1,help;
  int s_qu;
  
  //squark - Photon - squark
  Flavour flph = Flavour(kf_photon);
  if (flph.IsOn()) {
    //sUpypes
    for (short int i=1 ;i<7;i++) {
      if (i<4) s_qu=1000000 + 2*i;
      else     s_qu=2000000 + 2*i - 6;
      Flavour flav = Flavour((kf_code)(s_qu));
      if (flav.IsOn()) {
	vertex[vanz].in[0] = flav;
	vertex[vanz].in[1] = flph;
	vertex[vanz].in[2] = flav;
	
	Kabbala charge = Kabbala(string("Q_{"+flav.TexName()+"}"),flav.Charge());

	kcpl0 = -M_I*charge*g1;
	kcpl1 = kcpl0;

	
	vertex[vanz].cpl[0]  = kcpl0;
	vertex[vanz].cpl[1]  = kcpl1;
	vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	
	
	vertex[vanz].Color.push_back(Color_Function(cf::D));     
	vertex[vanz].Color.back().SetParticleArg(0,2);     
	vertex[vanz].Color.back().SetStringArg('0','2');     
	  
	
	vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("SSV",LF_Key()));
	vertex[vanz].Lorentz.back()->SetParticleArg(0,2,1);     
	
	vertex[vanz].on      = 1;
	vertex.push_back(Single_Vertex());vanz++;
      }  
    }

    //sDowntypes
    for (short int i=1 ;i<7;i++) {
      if (i<4) s_qu=1000000 + 2*i - 1;
      else     s_qu=2000000 + 2*i - 7;
      Flavour flav = Flavour((kf_code)(s_qu));
      if (flav.IsOn()) {
	vertex[vanz].in[0] = flav;
	vertex[vanz].in[1] = flph;
	vertex[vanz].in[2] = flav;
	
	Kabbala charge = Kabbala(string("Q_{"+flav.TexName()+"}"),flav.Charge());

	kcpl0 = -M_I*charge*g1;
	kcpl1 = kcpl0;
	
	vertex[vanz].cpl[0]  = kcpl0;
	vertex[vanz].cpl[1]  = kcpl1;
	vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	
	
	vertex[vanz].Color.push_back(Color_Function(cf::D));     
	vertex[vanz].Color.back().SetParticleArg(0,2);     
	vertex[vanz].Color.back().SetStringArg('0','2');     
	
	
	vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("SSV",LF_Key()));
	vertex[vanz].Lorentz.back()->SetParticleArg(0,2,1);     
	
	vertex[vanz].on      = 1;
	vertex.push_back(Single_Vertex());vanz++;
      }
    }
  }   

  //squark - Z - squark

  Flavour flZ = Flavour(kf_Z);
  if (flZ.IsOn()) {
    //sUptypes
    for (short int i=1;i<7;i++) {
      if (i<4) s_qu=1000000 + 2*i;
      else     s_qu=2000000 + 2*i - 6;
      Flavour flav1 = Flavour((kf_code)(s_qu));
      for (short int j=i;j<7;j++) {
	if (j<4) s_qu=1000000 + 2*j;
	else     s_qu=2000000 + 2*j - 6;
	Flavour flav2 = Flavour((kf_code)(s_qu));
	if (flav1.IsOn() && flav2.IsOn()) {
	  
	  vertex[vanz].in[0] = flav1;
	  vertex[vanz].in[1] = flZ;
	  vertex[vanz].in[2] = flav2;
	  
	  help = K_zero;
	  
	  if (i==j) {help = sintW*sintW*num_2/num_3;}  

	  kcpl0 = -M_I*g2/costW*
	    (K_Z_U(gen_sUp(flav2),j-1)*K_Z_U(gen_sUp(flav2),i-1)/num_2-help);
	  kcpl1 = kcpl0;
	  
	  vertex[vanz].cpl[0]  = kcpl0;
	  vertex[vanz].cpl[1]  = kcpl1;
	  vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	
	  
	  vertex[vanz].Color.push_back(Color_Function(cf::D));     
	  vertex[vanz].Color.back().SetParticleArg(0,2);     
	  vertex[vanz].Color.back().SetStringArg('0','2');     
	  
	  
	  vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("SSV",LF_Key()));
	  vertex[vanz].Lorentz.back()->SetParticleArg(0,2,1);     
	  
	  vertex[vanz].on      = 1;
	  vertex.push_back(Single_Vertex());vanz++;
	} 
      }
    }
    //sDowntypes
    for (short int i=1 ;i<7;i++) {
      if (i<4) s_qu=1000000 + 2*i - 1;
      else     s_qu=2000000 + 2*i - 7;
      Flavour flav1 = Flavour((kf_code)(s_qu));
      for (short int j=i ;j<7;j++) {
	if (j<4) s_qu=1000000 + 2*j - 1;
	else     s_qu=2000000 + 2*j - 7;
	Flavour flav2 = Flavour((kf_code)(s_qu));
	if (flav1.IsOn() && flav2.IsOn()) {
	  
	  vertex[vanz].in[0] = flav1;
	  vertex[vanz].in[1] = flZ;
	  vertex[vanz].in[2] = flav2;
	  
	  help = K_zero;
	  
	  if (i==j) { help = sintW*sintW/num_3; }  
	  
	  kcpl0 = M_I*g2/costW*
	    (K_Z_D(gen_sDown(flav2),j-1)*K_Z_D(gen_sDown(flav2),i-1)/num_2-help);
	  kcpl1 = kcpl0;

	  vertex[vanz].cpl[0]  = kcpl0;
	  vertex[vanz].cpl[1]  = kcpl1;
	  vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	  
	  
	  vertex[vanz].Color.push_back(Color_Function(cf::D));     
	  vertex[vanz].Color.back().SetParticleArg(0,2);     
	  vertex[vanz].Color.back().SetStringArg('0','2');     
	  
	  
	  vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("SSV",LF_Key()));
	  vertex[vanz].Lorentz.back()->SetParticleArg(0,2,1);     
	  
	  vertex[vanz].on      = 1;
	  vertex.push_back(Single_Vertex());vanz++;
	}  
      }
    }
  }    
  
  //supquarks - W - sdownquarks
  Flavour flW = Flavour(kf_Wplus);
  if (flW.IsOn()) {
    for (short int i=1;i<7;i++) {
      if (i<4) s_qu=1000000 + 2*i;
      else     s_qu=2000000 + 2*i - 6;
      Flavour flav1 = Flavour((kf_code)(s_qu));
      for (short int j=1;j<7;j++) {
	if (j<4) s_qu=1000000 + 2*j - 1;
	else     s_qu=2000000 + 2*j - 7;
	Flavour flav2 = Flavour((kf_code)(s_qu));
	if (flav1.IsOn() && flav2.IsOn()) {
	  
	  vertex[vanz].in[0] = flav1;
	  vertex[vanz].in[1] = flW;
	  vertex[vanz].in[2] = flav2;
	
	  Kabbala factor = K_zero;
	  
	  for (int I=0;I<3;I++) {
	    for (int J=0;J<3;J++) {
	      factor += K_Z_D(I,j-1)*K_Z_U(J,i-1)*
				    conj_K_CKM(J,I);
	    }
	  }

	  kcpl0 = -M_I*g2*invroot2*factor;
	  kcpl1 = kcpl0;

	  vertex[vanz].cpl[0]  = kcpl0;
	  vertex[vanz].cpl[1]  = kcpl1;
	  vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	  
	  
	  vertex[vanz].Color.push_back(Color_Function(cf::D));     
	  vertex[vanz].Color.back().SetParticleArg(0,2);     
	  vertex[vanz].Color.back().SetStringArg('0','2');     
	  
	  
	  vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("SSV",LF_Key()));
	  vertex[vanz].Lorentz.back()->SetParticleArg(0,2,1);     
	  
	  vertex[vanz].on      = 1;
	  vertex.push_back(Single_Vertex());vanz++;
	} 
      }			
    }
  } 
}

void Interaction_Model_sQuark_EW::c_SSS(std::vector<Single_Vertex>& vertex,int& vanz)
{
  Kabbala kcpl0,kcpl1,help;
  int s_qu;
  //sQuarks - A0 - sQuarks
  
  Flavour flA0 = Flavour(kf_A0);
  if (flA0.IsOn()) {
    //uptypes
    for (short int i=1;i<7;i++) {
      if (i<4) s_qu=1000000 + 2*i;
      else     s_qu=2000000 + 2*i - 6;
      Flavour flav1 = Flavour((kf_code)(s_qu));
      for (short int j=i;j<7;j++) {
	if (j<4) s_qu=1000000 + 2*j;
	else     s_qu=2000000 + 2*j - 6;
     	Flavour flav2 = Flavour((kf_code)(s_qu));
	if (flav1.IsOn() && flav2.IsOn()) {
	  
	  Kabbala K_uI = Kabbala(string("u^I"),Flavour((kf_code)(2*gen_sUp(flav1)+2)).Yuk()/
				 (v2).Value()*sqrt(2.));
	  
	  vertex[vanz].in[0] = flav1.Bar();
	  vertex[vanz].in[1] = flA0;
	  vertex[vanz].in[2] = flav2.Bar();
	  
	  kcpl0 = -(K_uI*K_Z_H(0,0)*(mu*K_Z_U(gen_sUp(flav1),i-1)*K_Z_U(gen_sUp(flav1)+3,j-1)-
				    conj_mu*K_Z_U(gen_sUp(flav1),j-1)*
				    K_Z_U(gen_sUp(flav1)+3,i-1))+
		   K_Z_H(1,0)*(K_u_S(gen_sUp(flav1),gen_sUp(flav2))*K_Z_U(gen_sUp(flav1),j-1)*
			       K_Z_U(gen_sUp(flav2)+3,i-1)-
			       K_u_S(gen_sUp(flav1),gen_sUp(flav2))*K_Z_U(gen_sUp(flav1),i-1)*
			       K_Z_U(gen_sUp(flav2)+3,j-1))+
		   K_Z_H(0,0)*(K_w_S(gen_sUp(flav1),gen_sUp(flav2))*K_Z_U(gen_sUp(flav1),i-1)*
			       K_Z_U(gen_sUp(flav2)+3,j-1)-K_w_S(gen_sUp(flav1),gen_sUp(flav2))*
			       K_Z_U(gen_sUp(flav1),j-1)*K_Z_U(gen_sUp(flav2)+3,i-1)))*invroot2;
	  
	  kcpl1 = kcpl0;
	  
	  vertex[vanz].cpl[0]  = kcpl0;
	  vertex[vanz].cpl[1]  = kcpl1;
	  vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	  
	  
	  vertex[vanz].Color.push_back(Color_Function(cf::D));     
	  vertex[vanz].Color.back().SetParticleArg(0,2);     
	  vertex[vanz].Color.back().SetStringArg('0','2');     
	  
	  
	  vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("SSS",LF_Key()));
	  
	  vertex[vanz].on      = 1;
	  vertex.push_back(Single_Vertex());vanz++;
	}
      }
    }
    //downtypes
    for (short int i=1;i<7;i++) {
      if (i<4) s_qu=1000000 + 2*i - 1;
      else     s_qu=2000000 + 2*i - 7;
      Flavour flav1 = Flavour((kf_code)(s_qu));
      for (short int j=i;j<7;j++) {
	if (j<4) s_qu=1000000 + 2*j - 1;
	else     s_qu=2000000 + 2*j - 7;
	Flavour flav2 = Flavour((kf_code)(s_qu));
	if (flav1.IsOn() && flav2.IsOn()) {
	  
	  Kabbala K_dI = Kabbala(string("d^I"),
				 -Flavour((kf_code)(2*gen_sDown(flav1)+1)).Yuk()/(v1).Value()*sqrt(2.));
	  
	  vertex[vanz].in[0] = flav1;
	  vertex[vanz].in[1] = flA0;
	  vertex[vanz].in[2] = flav2;
	  
	  kcpl0 = -(K_dI*K_Z_H(1,0)*(conj_mu*K_Z_D(gen_sDown(flav1),i-1)*
				     K_Z_D(gen_sDown(flav1)+3,j-1)
				     -mu*K_Z_D(gen_sDown(flav1),j-1)*
				     K_Z_D(gen_sDown(flav1)+3,i-1))
		   
		   + K_Z_H(0,0)*(K_d_S(gen_sDown(flav1),gen_sDown(flav2))*
				  K_Z_D(gen_sDown(flav1),j-1)*
				  K_Z_D(gen_sDown(flav2)+3,i-1)-
				  K_d_S(gen_sDown(flav1),gen_sDown(flav2))*
				  K_Z_D(gen_sDown(flav1),i-1)*
				  K_Z_D(gen_sDown(flav2)+3,j-1))
		   
		   +K_Z_H(1,0)*(K_e_S(gen_sDown(flav1),gen_sDown(flav2))*
				 K_Z_D(gen_sDown(flav1),j-1)*
				 K_Z_D(gen_sDown(flav2)+3,i-1)-
				 K_e_S(gen_sDown(flav1),gen_sDown(flav2))*
				 K_Z_D(gen_sDown(flav1),i-1)*
				 K_Z_D(gen_sDown(flav2)+3,j-1)))*invroot2;
	  
	  kcpl1 = kcpl0;
	  
	  vertex[vanz].cpl[0]  = kcpl0;
	  vertex[vanz].cpl[1]  = kcpl1;
	  vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	  
	  
	  vertex[vanz].Color.push_back(Color_Function(cf::D));     
	  vertex[vanz].Color.back().SetParticleArg(0,2);     
	  vertex[vanz].Color.back().SetStringArg('0','2');     
	  
	  
	  vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("SSS",LF_Key()));
	  
	  vertex[vanz].on      = 1;
	  vertex.push_back(Single_Vertex());vanz++;
	}
      }
    }
  } 
  
  
  //sQuarks - h0/H0 - sQuarks
  
  for (short int k=0;k<2;k++) {
    Flavour flH = Flavour((kf_code)(25+k*10)); 
    if (flH.IsOn()) {
    //uptypes  
    for (short int i=1;i<7;i++) {
      if (i<4) s_qu=1000000 + 2*i;
      else     s_qu=2000000 + 2*i - 6;
      Flavour flav1 = Flavour((kf_code)(s_qu));
      for (short int j=i;j<7;j++) {
	if (j<4) s_qu=1000000 + 2*j;
	else     s_qu=2000000 + 2*j - 6;
	Flavour flav2 =Flavour((kf_code)(s_qu));
	if(flav1.IsOn() && flav2.IsOn()){
	  vertex[vanz].in[0] = flav1;
	  vertex[vanz].in[1] = flH;
	  vertex[vanz].in[2] = flav2;
	  
	  help = K_zero;
	  
	  Kabbala K_uI = Kabbala(string("u^I"),Flavour((kf_code)(2*gen_sUp(flav1)+2)).Yuk()/
				 (v2).Value()*sqrt(2.));
	  
	  Kabbala fac = Kabbala(string("\\frac{3-8sin^2\\theta_W}{4sin^2\\theta_W}"),
				(3.-8.*(sintW).Value()*(sintW).Value())/
				(4.*(sintW).Value()*(sintW).Value()));
	  
	  if (i==j) {help = num_1;}
	  
	  kcpl0 = M_I*(-g1*g1/(costW*costW*num_3)*K_B_R(k)*
		       (help+fac*K_Z_U(gen_sUp(flav1),i-1)*K_Z_U(gen_sUp(flav1),j-1))
		       -K_uI*K_uI*v2*K_Z_R(1,k)*
		       (K_Z_U(gen_sUp(flav1),i-1)*K_Z_U(gen_sUp(flav1),j-1)+
		 	K_Z_U(gen_sUp(flav1)+3,i-1)*K_Z_U(gen_sUp(flav1)+3,j-1))
		       +K_Z_R(1,k)*invroot2*K_u_S(gen_sUp(flav1),gen_sUp(flav2))*
		       (K_Z_U(gen_sUp(flav1),i-1)*K_Z_U(gen_sUp(flav2)+3,j-1)+
			K_Z_U(gen_sUp(flav1),j-1)*K_Z_U(gen_sUp(flav2)+3,i-1))
		       +K_Z_R(0,k)*invroot2*K_w_S(gen_sUp(flav1),gen_sUp(flav2))*
		       (K_Z_U(gen_sUp(flav1),i-1)*K_Z_U(gen_sUp(flav2)+3,j-1)+
			K_Z_U(gen_sUp(flav1),j-1)*K_Z_U(gen_sUp(flav2)+3,i-1))
		       +K_uI*K_Z_R(0,k)*invroot2*mu*(K_Z_U(gen_sUp(flav1),j-1)*
							K_Z_U(gen_sUp(flav1)+3,i-1)+
							K_Z_U(gen_sUp(flav1),i-1)*
							K_Z_U(gen_sUp(flav1)+3,j-1))
		       );
	  
	  kcpl1 = kcpl0;
	  
	  vertex[vanz].cpl[0]  = kcpl0;
	  vertex[vanz].cpl[1]  = kcpl1;
	  vertex[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();
	  
	  
	  vertex[vanz].Color.push_back(Color_Function(cf::D));     
	  vertex[vanz].Color.back().SetParticleArg(0,2);     
	  vertex[vanz].Color.back().SetStringArg('0','2');     
	  
	  
	  vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("SSS",LF_Key()));
	  
	  vertex[vanz].on      = 1;
	  vertex.push_back(Single_Vertex());vanz++;
  	}
      }
    }
    
    //downtypes
    for (short int i=1;i<7;i++) {
      if (i<4) s_qu=1000000 + 2*i - 1;
      else     s_qu=2000000 + 2*i - 7;
      Flavour flav1 = Flavour((kf_code)(s_qu));
      for (short int j=i;j<7;j++) {
	if (j<4) s_qu=1000000 + 2*j - 1;
	else     s_qu=2000000 + 2*j - 7;
	Flavour flav2 =Flavour((kf_code)(s_qu));
	if(flav1.IsOn() && flav2.IsOn()){
	  vertex[vanz].in[0] = flav1;
	  vertex[vanz].in[1] = flH;
	  vertex[vanz].in[2] = flav2;
	  
	  help = K_zero;
	  
	  Kabbala K_dI = Kabbala(string("d^I"),
				 -Flavour((kf_code)(2*gen_sDown(flav1)+1)).Yuk()/(v1).Value()*sqrt(2.));
	  
	  Kabbala fac = Kabbala(string("\\frac{3-4sin^2\\theta_W}{2sin^2\\theta_W}"),
				(3.-4.*(sintW).Value()*(sintW).Value())/
				(2.*(sintW).Value()*(sintW).Value()));
	  
	  if (i==j) {help = num_1;}
	  
	  kcpl0 = M_I*(g1*g1/(costW*costW*num_6)*K_B_R(k)*
		       (help+fac*K_Z_D(gen_sDown(flav1),i-1)*K_Z_D(gen_sDown(flav1),j-1))
		       -K_dI*K_dI*v1*K_Z_R(0,k)*
		       (K_Z_D(gen_sDown(flav1),i-1)*K_Z_D(gen_sDown(flav1),j-1)+
			K_Z_D(gen_sDown(flav1)+3,i-1)*K_Z_D(gen_sDown(flav1)+3,j-1))
		       -K_Z_R(0,k)*invroot2*K_d_S(gen_sDown(flav1),gen_sDown(flav2))*
		       (K_Z_D(gen_sDown(flav1),j-1)*K_Z_D(gen_sDown(flav2)+3,i-1)+
			K_Z_D(gen_sDown(flav1),i-1)*K_Z_D(gen_sDown(flav2)+3,j-1))
		       +K_Z_R(1,k)*invroot2*K_e_S(gen_sDown(flav1),gen_sDown(flav2))*
		       (K_Z_D(gen_sDown(flav1),j-1)*K_Z_D(gen_sDown(flav2)+3,i-1)+
			K_Z_D(gen_sDown(flav1),i-1)*K_Z_D(gen_sDown(flav2)+3,j-1))
		       -K_dI*K_Z_R(1,k)*invroot2*(conj_mu*K_Z_D(gen_sDown(flav1),i-1)*
						     K_Z_D(gen_sDown(flav1)+3,j-1)+
						     mu*K_Z_D(gen_sDown(flav1),j-1)*
						     K_Z_D(gen_sDown(flav1)+3,i-1)));
	  
	  kcpl1 = kcpl0;
	  
	  vertex[vanz].cpl[0]  = kcpl0;
	  vertex[vanz].cpl[1]  = kcpl1;
	  vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	  
	  
	  vertex[vanz].Color.push_back(Color_Function(cf::D));     
	  vertex[vanz].Color.back().SetParticleArg(0,2);     
	  vertex[vanz].Color.back().SetStringArg('0','2');     
	  
	  
	  vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("SSS",LF_Key()));
	  
	  vertex[vanz].on      = 1;
	  vertex.push_back(Single_Vertex());vanz++;
	  
	}
      }
    }
    }
  }  
  //sUp - H+ - sDown
  
  Flavour flHplus = Flavour(kf_Hplus);
  if (flHplus.IsOn()) {
    for (short int i=1;i<7;i++) {
      if (i<4) s_qu=1000000 + 2*i;
      else     s_qu=2000000 + 2*i - 6;
      Flavour flav1 = Flavour((kf_code)(s_qu));
      for (short int j=1;j<7;j++) {
	if (j<4) s_qu=1000000 + 2*j - 1;
	else     s_qu=2000000 + 2*j - 7;
	Flavour flav2 = Flavour((kf_code)(s_qu));
	vertex[vanz].in[0] = flav1;
	vertex[vanz].in[1] = flHplus;
	vertex[vanz].in[2] = flav2;
	
	Kabbala K_dI = Kabbala(string("d^I"),
			       -Flavour((kf_code)(2*gen_sUp(flav1)+1)).Yuk()/(v1).Value()*sqrt(2.));
	
	Kabbala K_uJ = Kabbala(string("u^I"),Flavour((kf_code)(2*gen_sDown(flav2)+2)).Yuk()/
			       (v2).Value()*sqrt(2.));
	
	Kabbala K_massW = Kabbala(string("M_W"),Flavour(kf_Wplus).Mass());
	
	kcpl0 = M_I*((-(g2*g2)/num_2*(v1*K_Z_H(0,0)+v2*K_Z_H(1,0))+
		      v1*K_dI*K_dI*K_Z_H(0,0)+v2*K_uJ*K_uJ*K_Z_H(1,0))*invroot2*
		     K_CKM(gen_sUp(flav1),gen_sDown(flav2))*
		     K_Z_D(gen_sDown(flav2),j-1)*K_Z_U(gen_sUp(flav1),i-1)

		     - sintW*K_massW*root2/g1*K_uJ*K_dI*
		     K_CKM(gen_sUp(flav1),gen_sDown(flav2))*
		     K_Z_D(gen_sDown(flav2)+3,j-1)*K_Z_U(gen_sUp(flav1)+3,i-1)

		     + (K_Z_H(0,0)*conj_mu*K_uJ*K_CKM(gen_sUp(flav1),gen_sDown(flav2))
			+(K_Z_H(0,0)*K_w_S(0,gen_sUp(flav1))-K_Z_H(1,0)*K_u_S(0,gen_sUp(flav1)))*
			conj_K_CKM(gen_sUp(flav1),0)+
			(K_Z_H(0,0)*K_w_S(1,gen_sUp(flav1))-K_Z_H(1,0)*K_u_S(1,gen_sUp(flav1)))*
			conj_K_CKM(gen_sUp(flav1),1)+
			(K_Z_H(0,0)*K_w_S(2,gen_sUp(flav1))-K_Z_H(1,0)*K_u_S(2,gen_sUp(flav1)))*
			K_CKM(gen_sUp(flav1),2))*
		     K_Z_U(gen_sUp(flav1)+3,i-1)*K_Z_D(gen_sDown(flav2),j-1)
		     + ((K_Z_H(0,0)*K_d_S(0,gen_sUp(flav1))+K_Z_H(1,0)*K_e_S(0,gen_sUp(flav1)))*
			conj_K_CKM(0,gen_sDown(flav2))+
			(K_Z_H(0,0)*K_d_S(1,gen_sUp(flav1))+K_Z_H(1,0)*K_e_S(1,gen_sUp(flav1)))*
			conj_K_CKM(1,gen_sDown(flav2))+
			(K_Z_H(0,0)*K_d_S(2,gen_sUp(flav1))+K_Z_H(1,0)*K_e_S(2,gen_sUp(flav1)))*
			conj_K_CKM(2,gen_sDown(flav2))-
			K_Z_H(1,0)*mu*K_dI*conj_K_CKM(gen_sUp(flav1),gen_sDown(flav2)))*
		     K_Z_U(gen_sDown(flav2),i-1)*K_Z_D(gen_sUp(flav1)+3,j-1));
	
	kcpl1 = kcpl0;
	
	vertex[vanz].cpl[0]  = kcpl0;
	vertex[vanz].cpl[1]  = kcpl1;
	vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	
	
	vertex[vanz].Color.push_back(Color_Function(cf::D));     
	vertex[vanz].Color.back().SetParticleArg(0,2);     
	vertex[vanz].Color.back().SetStringArg('0','2');     
	
	
	vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("SSS",LF_Key()));
	
	vertex[vanz].on      = 1;
	vertex.push_back(Single_Vertex());vanz++;
      }
    }
  }
}

void Interaction_Model_sQuark_EW::c_SSVV(std::vector<Single_Vertex>& vertex,int& vanz)
{
  Flavour flavW(kf_Wplus);
  Flavour flavZ(kf_Z);
  Flavour flavPhoton(kf_photon);
  Kabbala kcpl0,kcpl1,help;
  int s_qu;
  
  for (short int l=1;l<3;l++) {
    for (short int i=1;i<7;i++) {
      s_qu = l*1000000 + i;
      Flavour flav = Flavour((kf_code)(s_qu));
      if (flav.IsOn() && flavPhoton.IsOn()) {
	// P - U/D - U/D - P  
	vertex[vanz].in[0] = flavPhoton;
	vertex[vanz].in[1] = flav.Bar();
	vertex[vanz].in[2] = flav;
	vertex[vanz].in[3] = flavPhoton;
	
	vertex[vanz].nleg     = 4;
	
	 Kabbala charge = Kabbala(string("Q_{"+flav.TexName()+"}"),flav.Charge());
	 
	 kcpl0 = M_I*charge*charge*num_2*g1*g1;
	 kcpl1 = kcpl0;
	 
	 vertex[vanz].cpl[0]  = kcpl0;
	 vertex[vanz].cpl[1]  = kcpl1;
	 vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	 
	 
	 vertex[vanz].Color.push_back(Color_Function(cf::D));     
	 vertex[vanz].Color.back().SetParticleArg(1,2);     
	 vertex[vanz].Color.back().SetStringArg('1','2');     
	 
	 
	 vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("VVSS",LF_Key()));
	 vertex[vanz].Lorentz.back()->SetParticleArg(0,3);     
	 
	 vertex[vanz].on      = 1;
	 vertex.push_back(Single_Vertex());vanz++;
      }
    }
  }
  
  // P - U - U - Z  
  for (short int i=1;i<7;i++) {
    if (i<4) s_qu=1000000 + 2*i;
    else     s_qu=2000000 + 2*i - 6;
    Flavour flav1 = Flavour((kf_code)(s_qu));
    for (short int j=1;j<7;j++) {
      if (j<4) s_qu=1000000 + 2*j;
      else     s_qu=2000000 + 2*j - 6;
      Flavour flav2 = Flavour((kf_code)(s_qu));
      if (flavPhoton.IsOn() && flavZ.IsOn()) {
	if (flav1.IsOn() && flav2.IsOn() && gen_sUp(flav1)==gen_sUp(flav2)) {

	  vertex[vanz].in[0] = flavPhoton;
	  vertex[vanz].in[1] = flav1.Bar();
	  vertex[vanz].in[2] = flav2;
	  vertex[vanz].in[3] = flavZ;
	  
	  vertex[vanz].nleg     = 4;

	  help = K_zero;
	  
	  if (i==j) {help = sintW*sintW*num_4/num_3;}  

	  kcpl0 = num_2*M_I*g1*g2/(num_3*costW)*
	    ((K_Z_U(0,j-1)*K_Z_U(0,i-1)+ 
	      K_Z_U(1,j-1)*K_Z_U(1,i-1)+ 
	      K_Z_U(2,j-1)*K_Z_U(2,i-1)) - help);
	  kcpl1 = kcpl0;
	  	  
	  vertex[vanz].cpl[0]  = kcpl0;
	  vertex[vanz].cpl[1]  = kcpl1;
	  vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	  
	  
	  vertex[vanz].Color.push_back(Color_Function(cf::D));     
	  vertex[vanz].Color.back().SetParticleArg(1,2);     
	  vertex[vanz].Color.back().SetStringArg('1','2');     
	  
	  
	  vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("VVSS",LF_Key()));
	  vertex[vanz].Lorentz.back()->SetParticleArg(0,3);     
	  
	  vertex[vanz].on      = 1;
	  vertex.push_back(Single_Vertex());vanz++;
	}
      }
    }
  }
  // P - D - D - Z  
  for (short int i=1;i<7;i++) {
    if (i<4) s_qu=1000000 + 2*i - 1;
    else     s_qu=2000000 + 2*i - 7;
    Flavour flav1 = Flavour((kf_code)(s_qu));
    for (short int j=1;j<7;j++) {
      if (j<4) s_qu=1000000 + 2*j - 1;
      else     s_qu=2000000 + 2*j - 7;
      Flavour flav2 = Flavour((kf_code)(s_qu));
      if (flavPhoton.IsOn() && flavZ.IsOn()) {
	if (flav1.IsOn() && flav2.IsOn() && gen_sDown(flav1)==gen_sDown(flav2)) {

	  vertex[vanz].in[0] = flavPhoton;
	  vertex[vanz].in[1] = flav1.Bar();
	  vertex[vanz].in[2] = flav2;
	  vertex[vanz].in[3] = flavZ;
	  
	  vertex[vanz].nleg     = 4;

	  help = K_zero;
	  
	  if (i==j) {help = sintW*sintW*num_2/num_3;}  

	  kcpl0 = M_I*g1*g2/(num_3*costW)*
	    ((K_Z_D(0,j-1)*K_Z_D(0,i-1) +
	      K_Z_D(1,j-1)*K_Z_D(1,i-1) +
	      K_Z_D(2,j-1)*K_Z_D(2,i-1)) - help);
	  kcpl1 = kcpl0;
	  	  
	  vertex[vanz].cpl[0]  = kcpl0;
	  vertex[vanz].cpl[1]  = kcpl1;
	  vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	  
	  
	  vertex[vanz].Color.push_back(Color_Function(cf::D));     
	  vertex[vanz].Color.back().SetParticleArg(1,2);     
	  vertex[vanz].Color.back().SetStringArg('1','2');     
	  
	  
	  vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("VVSS",LF_Key()));
	  vertex[vanz].Lorentz.back()->SetParticleArg(0,3);     
	  
	  vertex[vanz].on      = 1;
	  vertex.push_back(Single_Vertex());vanz++;
	}
      }
    }
  }

  // Z - U - U - Z  
  for (short int i=1;i<7;i++) {
    if (i<4) s_qu=1000000 + 2*i;
    else     s_qu=2000000 + 2*i - 6;
    Flavour flav1 = Flavour((kf_code)(s_qu));
    for (short int j=1;j<7;j++) {
      if (j<4) s_qu=1000000 + 2*j;
      else     s_qu=2000000 + 2*j - 6;
      Flavour flav2 = Flavour((kf_code)(s_qu));
      if (flavZ.IsOn()) {
	if (flav1.IsOn() && flav2.IsOn() && gen_sUp(flav1)==gen_sUp(flav2)) {
	  
	  help = K_zero;
	  
	  if (i==j) {help = num_1;}  
	  
	  vertex[vanz].in[0] = flavZ;
	  vertex[vanz].in[1] = flav1.Bar();
	  vertex[vanz].in[2] = flav2;
	  vertex[vanz].in[3] = flavZ;
	  
	  vertex[vanz].nleg     = 4;
	  
	  kcpl0 = num_2*M_I*g1*g1/(num_3*costW*costW)*
	    (num_4/num_3*help*sintW*sintW + (num_3 - num_2*num_4*sintW*sintW)/(num_4*sintW*sintW)*
	     (K_Z_U(0,i-1)*K_Z_U(0,j-1) + 
	      K_Z_U(1,i-1)*K_Z_U(1,j-1) + 
	      K_Z_U(2,i-1)*K_Z_U(2,j-1)));
	  kcpl1 = kcpl0;
	  
	  vertex[vanz].cpl[0]  = kcpl0;
	  vertex[vanz].cpl[1]  = kcpl1;
	  vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	  
	  
	  vertex[vanz].Color.push_back(Color_Function(cf::D));     
	  vertex[vanz].Color.back().SetParticleArg(1,2);     
	  vertex[vanz].Color.back().SetStringArg('1','2');     
	  
	  
	  vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("VVSS",LF_Key()));
	  vertex[vanz].Lorentz.back()->SetParticleArg(0,3);     
	  
	  vertex[vanz].on      = 1;
	  vertex.push_back(Single_Vertex());vanz++;
	}
      }
    }
  }
  // Z - D - D - Z  
  for (short int i=1;i<7;i++) {
    if (i<4) s_qu=1000000 + 2*i - 1;
    else     s_qu=2000000 + 2*i - 7;
    Flavour flav1 = Flavour((kf_code)(s_qu));
    for (short int j=1;j<7;j++) {
      if (j<4) s_qu=1000000 + 2*j - 1;
      else     s_qu=2000000 + 2*j - 7;
      Flavour flav2 = Flavour((kf_code)(s_qu));
      if (flavZ.IsOn()) {
	if (flav1.IsOn() && flav2.IsOn() && gen_sDown(flav1)==gen_sDown(flav2)) {
	  
	  help = K_zero;
	  
	  if (i==j) {help = num_1;}  
	  
	  vertex[vanz].in[0] = flavZ;
	  vertex[vanz].in[1] = flav1.Bar();
	  vertex[vanz].in[2] = flav2;
	  vertex[vanz].in[3] = flavZ;
	  
	  vertex[vanz].nleg     = 4;
	  
	  kcpl0 = num_2*M_I*g1*g1/(num_3*costW*costW)*
	    (num_1/num_3*help*sintW*sintW + (num_3 - num_4*sintW*sintW)/(num_4*sintW*sintW)*
	     (K_Z_D(0,i-1)*K_Z_D(0,j-1) + 
	      K_Z_D(1,i-1)*K_Z_D(1,j-1) + 
	      K_Z_D(2,i-1)*K_Z_D(2,j-1)));
	  kcpl1 = kcpl0;
	  
	  vertex[vanz].cpl[0]  = kcpl0;
	  vertex[vanz].cpl[1]  = kcpl1;
	  vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	  
	  
	  vertex[vanz].Color.push_back(Color_Function(cf::D));     
	  vertex[vanz].Color.back().SetParticleArg(1,2);     
	  vertex[vanz].Color.back().SetStringArg('1','2');     
	  
	  
	  vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("VVSS",LF_Key()));
	  vertex[vanz].Lorentz.back()->SetParticleArg(0,3);     
	  
	  vertex[vanz].on      = 1;
	  vertex.push_back(Single_Vertex());vanz++;
	}
      }
    }      
  }

  // W - D - U - P/Z  
  for (short int i=1;i<7;i++) {
    if (i<4) s_qu=1000000 + 2*i;
    else     s_qu=2000000 + 2*i - 6;
    Flavour flav1 = Flavour((kf_code)(s_qu));
    for (short int j=1;j<7;j++) {
      if (j<4) s_qu=1000000 + 2*j - 1;
      else     s_qu=2000000 + 2*j - 7;
      Flavour flav2 = Flavour((kf_code)(s_qu));
      if (flav1.IsOn() && flav2.IsOn()) {
	// W - D - U - P  
	if (flavW.IsOn()) {
	  if (flavPhoton.IsOn()) {

	    vertex[vanz].in[0] = flavW.Bar();
	    vertex[vanz].in[1] = flav1.Bar();
	    vertex[vanz].in[2] = flav2;
	    vertex[vanz].in[3] = flavPhoton;
	    
	    vertex[vanz].nleg     = 4;

	    Kabbala factor = K_zero;
	    
	    for (int I=0;I<3;I++) {
	      for (int J=0;J<3;J++) {
		factor += K_Z_D(I,j-1)*K_Z_U(J,i-1)*
				      conj_K_CKM(J,I);
	      }
	    }
	    
	    kcpl0 = M_I*g1*g2*root2/num_6*factor;
					   
	    kcpl1 = kcpl0;
	    
	    vertex[vanz].cpl[0]  = kcpl0;
	    vertex[vanz].cpl[1]  = kcpl1;
	    vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	    
	    
	    vertex[vanz].Color.push_back(Color_Function(cf::D));     
	    vertex[vanz].Color.back().SetParticleArg(1,2);     
	    vertex[vanz].Color.back().SetStringArg('1','2');     
	    
	    
	    vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("VVSS",LF_Key()));
	    vertex[vanz].Lorentz.back()->SetParticleArg(0,3);     
	    
	    vertex[vanz].on      = 1;
	    vertex.push_back(Single_Vertex());vanz++;
	  }
	  
	  if (flavZ.IsOn()) {
	    // W - D - U - Z  
	    
	    vertex[vanz].in[0] = flavW.Bar();
	    vertex[vanz].in[1] = flav1.Bar();
	    vertex[vanz].in[2] = flav2;
	    vertex[vanz].in[3] = flavZ;
	    
	    vertex[vanz].nleg     = 4;
	
	    Kabbala factor = K_zero;
	    
	    for (int I=0;I<3;I++) {
	      for (int J=0;J<3;J++) {
		factor += K_Z_D(I,j-1)*K_Z_U(J,i-1)*
				      conj_K_CKM(J,I);
	      }
	    }
	    
	    kcpl0 = -M_I*g1*g1*root2/(num_6*costW)*factor;
     	    
	    kcpl1 = kcpl0;
	    
	    vertex[vanz].cpl[0]  = kcpl0;
	    vertex[vanz].cpl[1]  = kcpl1;
	    vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	    
	    
	    vertex[vanz].Color.push_back(Color_Function(cf::D));     
	    vertex[vanz].Color.back().SetParticleArg(1,2);     
	    vertex[vanz].Color.back().SetStringArg('1','2');     
	    
	    
	    vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("VVSS",LF_Key()));
	    vertex[vanz].Lorentz.back()->SetParticleArg(0,3);     
	    
	    vertex[vanz].on      = 1;
	    vertex.push_back(Single_Vertex());vanz++;
	  }
	}
      }
    }
  }
  // W - D/U - D/U - W  
  if (flavW.IsOn()) {
    for (short int i=1;i<7;i++) {
      if (i<4) s_qu=1000000 + 2*i - 1;
      else     s_qu=2000000 + 2*i - 7;
      Flavour flav1 = Flavour((kf_code)(s_qu));
      for (short int j=1;j<7;j++) {
	if (j<4) s_qu=1000000 + 2*j - 1;
	else     s_qu=2000000 + 2*j - 7;
	Flavour flav2 = Flavour((kf_code)(s_qu));
	if (flav1.IsOn() && flav2.IsOn()) {
	  // W - D - D - W  
	  vertex[vanz].in[0] = flavW.Bar();
	  vertex[vanz].in[1] = flav1.Bar();
	  vertex[vanz].in[2] = flav2;
	  vertex[vanz].in[3] = flavW.Bar();
	    
	  vertex[vanz].nleg     = 4;
	  
	  kcpl0 = M_I*g2*g2/num_2*(K_Z_D(0,i-1)*K_Z_D(0,j-1) +
				   K_Z_D(1,i-1)*K_Z_D(1,j-1) +
				   K_Z_D(2,i-1)*K_Z_D(2,j-1));
	  kcpl1 = kcpl0;
	  
	  vertex[vanz].cpl[0]  = kcpl0;
	  vertex[vanz].cpl[1]  = kcpl1;
	  vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	  
	  
	  vertex[vanz].Color.push_back(Color_Function(cf::D));     
	  vertex[vanz].Color.back().SetParticleArg(1,2);     
	  vertex[vanz].Color.back().SetStringArg('1','2');     
	  
	  
	  vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("VVSS",LF_Key()));
	  vertex[vanz].Lorentz.back()->SetParticleArg(0,3);     
	  
	  vertex[vanz].on      = 1;
	  vertex.push_back(Single_Vertex());vanz++;
	}
      }
    }
    for (short int i=1;i<7;i++) {
      if (i<4) s_qu=1000000 + 2*i;
      else     s_qu=2000000 + 2*i - 6;
      Flavour flav1 = Flavour((kf_code)(s_qu));
      for (short int j=1;j<7;j++) {
	if (j<4) s_qu=1000000 + 2*j;
	else     s_qu=2000000 + 2*j - 6;
	Flavour flav2 = Flavour((kf_code)(s_qu));
	if (flav1.IsOn() && flav2.IsOn()) {
	  // W - U - U - W  
	  
	  vertex[vanz].in[0] = flavW.Bar();
	  vertex[vanz].in[1] = flav1.Bar();
	  vertex[vanz].in[2] = flav2;
	  vertex[vanz].in[3] = flavW.Bar();
	    
	  vertex[vanz].nleg     = 4;
	  
	  kcpl0 = M_I*g2*g2/num_2*(K_Z_U(0,i-1)*K_Z_U(0,j-1) +
				   K_Z_U(1,i-1)*K_Z_U(1,j-1) +
				   K_Z_U(2,i-1)*K_Z_U(2,j-1));

	  kcpl1 = kcpl0;
	  
	  vertex[vanz].cpl[0]  = kcpl0;
	  vertex[vanz].cpl[1]  = kcpl1;
	  vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	  
	  
	  vertex[vanz].Color.push_back(Color_Function(cf::D));     
	  vertex[vanz].Color.back().SetParticleArg(1,2);     
	  vertex[vanz].Color.back().SetStringArg('1','2');     
	  
	  
	  vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("VVSS",LF_Key()));
	  vertex[vanz].Lorentz.back()->SetParticleArg(0,3);     
	  
	  vertex[vanz].on      = 1;
	  vertex.push_back(Single_Vertex());vanz++;
	}
      }
    }	
  }

  //Interactions of two squarks an EW gauge boson and a gluon
  
  Flavour flgluon = Flavour(kf_gluon);
  
  if (flgluon.IsOn()) {
    // G - U/D - U/D - P  
    for (short int i=1;i<7;i++) {
      if (i<4) s_qu=1000000 + 2*i;
      else     s_qu=2000000 + 2*i - 6;
      Flavour flav = Flavour((kf_code)(s_qu));
      if (flav.IsOn() && flavPhoton.IsOn()) {
	// G - U - U - P  
 	vertex[vanz].in[0] = flgluon;
	vertex[vanz].in[1] = flav.Bar();
	vertex[vanz].in[2] = flav;
	vertex[vanz].in[3] = flavPhoton;
	
	vertex[vanz].nleg     = 4;
	
	Kabbala charge = Kabbala(string("Q_{"+flav.TexName()+"}"),flav.Charge());
	
	kcpl0 = M_I*charge*num_2*g1*g3;
	kcpl1 = kcpl0;
	
	vertex[vanz].cpl[0]  = kcpl0;
	vertex[vanz].cpl[1]  = kcpl1;
	vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	
	
	vertex[vanz].Color.push_back(Color_Function(cf::T));     
	vertex[vanz].Color.back().SetParticleArg(0,2,1);     
	vertex[vanz].Color.back().SetStringArg('0','2','1');    
	
	
	vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("VVSS",LF_Key()));
	vertex[vanz].Lorentz.back()->SetParticleArg(0,3);     
	
	vertex[vanz].on      = 1;
	vertex.push_back(Single_Vertex());vanz++;
      }
    }
    for (short int i=1;i<7;i++) {
      if (i<4) s_qu=1000000 + 2*i - 1;
      else     s_qu=2000000 + 2*i - 7;
      Flavour flav = Flavour((kf_code)(s_qu));
      if (flav.IsOn() && flavPhoton.IsOn()) {
	// G - D - D - P  
	vertex[vanz].in[0] = flgluon;
	vertex[vanz].in[1] = flav.Bar();
	vertex[vanz].in[2] = flav;
	vertex[vanz].in[3] = flavPhoton;
	
	vertex[vanz].nleg     = 4;
	
	Kabbala charge = Kabbala(string("Q_{"+flav.TexName()+"}"),flav.Charge());
	
	kcpl0 = M_I*charge*num_2*g1*g3;
	kcpl1 = kcpl0;
	
	vertex[vanz].cpl[0]  = kcpl0;
	vertex[vanz].cpl[1]  = kcpl1;
	vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	
	
	vertex[vanz].Color.push_back(Color_Function(cf::T));     
	vertex[vanz].Color.back().SetParticleArg(0,2,1);     
	vertex[vanz].Color.back().SetStringArg('0','2','1');    
		
	
	vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("VVSS",LF_Key()));
	vertex[vanz].Lorentz.back()->SetParticleArg(0,3);     
	
	vertex[vanz].on      = 1;
	vertex.push_back(Single_Vertex());vanz++;
      }
    }
    // G - U/D - U/D - Z  
    for (short int i=1;i<7;i++) {
      if (i<4) s_qu=1000000 + 2*i;
      else     s_qu=2000000 + 2*i - 6;
      Flavour flav1 = Flavour((kf_code)(s_qu));
      for (short int j=1;j<7;j++) {
	if (j<4) s_qu=1000000 + 2*j;
	else     s_qu=2000000 + 2*j - 6;
	Flavour flav2 = Flavour((kf_code)(s_qu));
	if (flav1.IsOn() && flav2.IsOn() && flavZ.IsOn()) {
	  if (gen_sUp(flav1)==gen_sUp(flav2)) {
	    // G - U - U - Z  
	    vertex[vanz].in[0] = flgluon;
	    vertex[vanz].in[1] = flav1.Bar();
	    vertex[vanz].in[2] = flav2;
	    vertex[vanz].in[3] = flavZ;
	    
	    vertex[vanz].nleg     = 4;
	    
	    help = K_zero;
	    
	    if (i==j) {help = sintW*sintW*num_4/num_3;}  
	    
	    kcpl0 = M_I*g2*g3/costW*
	      (K_Z_U(gen_sUp(flav2),j-1)*K_Z_U(gen_sUp(flav1),i-1) - help);
	    kcpl1 = kcpl0;
	    
	    vertex[vanz].cpl[0]  = kcpl0;
	    vertex[vanz].cpl[1]  = kcpl1;
	    vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	    
	    
	    vertex[vanz].Color.push_back(Color_Function(cf::T));     
	    vertex[vanz].Color.back().SetParticleArg(0,2,1);     
	    vertex[vanz].Color.back().SetStringArg('0','2','1');    
	    
	    
	    vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("VVSS",LF_Key()));
	    vertex[vanz].Lorentz.back()->SetParticleArg(0,3);     
	    
	    vertex[vanz].on      = 1;
	    vertex.push_back(Single_Vertex());vanz++;
	  }
	}
      }
    }
    for (short int i=1;i<7;i++) {
      if (i<4) s_qu=1000000 + 2*i - 1;
      else     s_qu=2000000 + 2*i - 7;
      Flavour flav1 = Flavour((kf_code)(s_qu));
      for (short int j=1;j<7;j++) {
	if (j<4) s_qu=1000000 + 2*j - 1;
	else     s_qu=2000000 + 2*j - 7;
	Flavour flav2 = Flavour((kf_code)(s_qu));
	if (flav1.IsOn() && flav2.IsOn() && flavZ.IsOn()) {
	  if (gen_sDown(flav1)==gen_sDown(flav2)) {
	    // G - D - D - Z  
	    vertex[vanz].in[0] = flgluon;
	    vertex[vanz].in[1] = flav1.Bar();
	    vertex[vanz].in[2] = flav2;
	    vertex[vanz].in[3] = flavZ;
	    
	    vertex[vanz].nleg     = 4;
	    
	    help = K_zero;
	    
	    if (i==j) {help = sintW*sintW*num_2/num_3;}  
	    
	    kcpl0 = M_I*g2*g3/costW*
	      (-K_Z_D(gen_sDown(flav1),j-1)*K_Z_D(gen_sDown(flav1),i-1) + help);
	    kcpl1 = kcpl0;
	    
	    vertex[vanz].cpl[0]  = kcpl0;
	    vertex[vanz].cpl[1]  = kcpl1;
	    vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	    
	    
	    vertex[vanz].Color.push_back(Color_Function(cf::T));     
	    vertex[vanz].Color.back().SetParticleArg(0,2,1);     
	    vertex[vanz].Color.back().SetStringArg('0','2','1');    
	    
	    
	    vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("VVSS",LF_Key()));
	    vertex[vanz].Lorentz.back()->SetParticleArg(0,3);     
	    
	    vertex[vanz].on      = 1;
	    vertex.push_back(Single_Vertex());vanz++;
	  }
	}
      }
    }

    // G - D - U - W  
    for (short int i=1;i<7;i++) {
      if (i<4) s_qu=1000000 + 2*i;
      else     s_qu=2000000 + 2*i - 6;
      Flavour flav1 = Flavour((kf_code)(s_qu));
      for (short int j=1;j<7;j++) {
	if (j<4) s_qu=1000000 + 2*j - 1;
	else     s_qu=2000000 + 2*j - 7;
	Flavour flav2 = Flavour((kf_code)(s_qu));
	if (flav1.IsOn() && flav2.IsOn()) {
	  
	  if (flavW.IsOn()) {
	    
	    vertex[vanz].in[0] = flavW.Bar();
	    vertex[vanz].in[1] = flav1.Bar();
	    vertex[vanz].in[2] = flav2;
	    vertex[vanz].in[3] = flgluon;
	    
	    vertex[vanz].nleg     = 4;
	    kcpl0 = M_I*g2*g3*root2*(K_Z_D(gen_sDown(flav2),j-1)*K_Z_U(gen_sUp(flav1),i-1))*
	      conj_K_CKM(gen_sUp(flav1),gen_sDown(flav2));
	    
	    kcpl1 = kcpl0;
	    
	    vertex[vanz].cpl[0]  = kcpl0;
	    vertex[vanz].cpl[1]  = kcpl1;
	    vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	    
	    
	    vertex[vanz].Color.push_back(Color_Function(cf::T));
	    vertex[vanz].Color.back().SetParticleArg(3,2,1);     
	    vertex[vanz].Color.back().SetStringArg('3','2','1');    
	    
	    
	    vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("VVSS",LF_Key()));
	    vertex[vanz].Lorentz.back()->SetParticleArg(0,3);     
	    
	    vertex[vanz].on      = 1;
	    vertex.push_back(Single_Vertex());vanz++;
	  }
	}
      }
    }
  }
}

void Interaction_Model_sQuark_EW::c_SSSS(std::vector<Single_Vertex>& vertex,int& vanz)
{
  Flavour flHmin=Flavour(kf_Hplus).Bar();    
  Flavour flA0(kf_A0);    
  Kabbala kcpl0,kcpl1,num_1,help;
  int s_qu;

  //up-squarks Higgs interactions
  for (short int i=1;i<7;++i) {
    if (i<4) s_qu=1000000 + 2*i;
    else     s_qu=2000000 + 2*i - 6;
    Flavour flav1 = Flavour((kf_code)(s_qu));
    for (short int j=i;j<7;++j) {
      if (j<4) s_qu=1000000 + 2*j;
      else     s_qu=2000000 + 2*j - 6;
      Flavour flav2 = Flavour((kf_code)(s_qu));
      if (flav1.IsOn() && flav2.IsOn() && gen_sUp(flav1)==gen_sUp(flav2)) {
	
	Kabbala K_uI = Kabbala(string("u^I"),Flavour((kf_code)(2*gen_sUp(flav1)+2)).Yuk()/
			       (v2).Value()*sqrt(2.));
	
	Kabbala K_dI = Kabbala(string("d^I"),
			       -Flavour((kf_code)(2*gen_sUp(flav1)+1)).Yuk()/(v1).Value()*sqrt(2.));
	
	//Hmin -> sup - sup - Hmin 
	if (flHmin.IsOn()) {
	  vertex[vanz].in[0] = flHmin;
	  vertex[vanz].in[1] = flav1.Bar();
	  vertex[vanz].in[2] = flav2;
	  vertex[vanz].in[3] = flHmin;
	  
	  vertex[vanz].nleg  = 4;  
	  
	  help = K_zero;
	  if(i==j) help = num_1;
	  
	  kcpl0 = M_I*(-g1*g1/(costW*costW*num_3)*K_A_H(0,0)*
		       (help-(num_3+num_2*sintW*sintW)/(num_4*sintW*sintW)*
			K_Z_U(gen_sUp(flav1),i-1)*K_Z_U(gen_sUp(flav1),j-1))
		       -K_uI*K_uI*K_Z_H(1,0)*K_Z_H(1,0)*
		       K_Z_U(gen_sUp(flav1)+3,i-1)*K_Z_U(gen_sUp(flav1)+3,j-1)
		       -K_dI*K_dI*K_Z_H(0,0)*K_Z_H(0,0)*
		       (K_Z_U(0,i-1)*K_CKM(0,gen_sUp(flav1)) +
			K_Z_U(1,i-1)*K_CKM(1,gen_sUp(flav1)) +
			K_Z_U(2,i-1)*K_CKM(1,gen_sUp(flav1)))*
		       (K_Z_U(0,j-1)*conj_K_CKM(0,gen_sUp(flav1)) +
			K_Z_U(1,j-1)*conj_K_CKM(1,gen_sUp(flav1)) +
			K_Z_U(2,j-1)*conj_K_CKM(1,gen_sUp(flav1))));
	  
	  kcpl1 = kcpl0;
	  
	  vertex[vanz].cpl[0]  = kcpl0; 
	  vertex[vanz].cpl[1]  = kcpl1;
	  vertex[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();
	  
	  
	  vertex[vanz].Color.push_back(Color_Function(cf::D));     
	  vertex[vanz].Color.back().SetParticleArg(1,2);     
	  vertex[vanz].Color.back().SetStringArg('1','2');     
	  
	  
	  vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("SSSS",LF_Key()));     
	  
	  vertex[vanz].on      = 1;
	  vertex.push_back(Single_Vertex());vanz++;
	}
	//A0 -> sUp - sUp - A0 
	if (flA0.IsOn()) {
	  vertex[vanz].in[0] = flA0;
	  vertex[vanz].in[1] = flav1.Bar();
	  vertex[vanz].in[2] = flav2;
	  vertex[vanz].in[3] = flA0;
	  
	  vertex[vanz].nleg  = 4;  
	  
	  help = K_zero;
	  if(i==j) help = num_1;
	  
	  
	  kcpl0 = M_I*(-g1*g1/(costW*costW*num_3)*K_A_H(0,0)*
		       (help+(num_3-num_4*num_2*sintW*sintW)/(num_4*sintW*sintW)*
			K_Z_U(gen_sUp(flav1),i-1)*K_Z_U(gen_sUp(flav1),j-1))
		       -K_uI*K_uI*K_Z_H(1,0)*K_Z_H(1,0)*
		       (K_Z_U(gen_sUp(flav1),i-1)*K_Z_U(gen_sUp(flav1),j-1) + 
			K_Z_U(gen_sUp(flav1)+3,i-1)*K_Z_U(gen_sUp(flav1)+3,j-1)));
	  
	  kcpl1 = kcpl0;
	  
	  vertex[vanz].cpl[0]  = kcpl0; 
	  vertex[vanz].cpl[1]  = kcpl1;
	  vertex[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();
	  
	  
	  vertex[vanz].Color.push_back(Color_Function(cf::D));     
	  vertex[vanz].Color.back().SetParticleArg(1,2);     
	  vertex[vanz].Color.back().SetStringArg('1','2');     
	    
	  
	  vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("SSSS",LF_Key()));     
	  
	  vertex[vanz].on      = 1;
	  vertex.push_back(Single_Vertex());vanz++;
	}
	//h0/H0 -> sUp - sUp - h0/H0 
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
	      
	      kcpl0 = M_I*(-g1*g1/(costW*costW*num_3)*K_A_R(k,l)*
			   (help+(num_3-num_4*num_2*sintW*sintW)/(num_4*sintW*sintW)*
			    K_Z_U(gen_sUp(flav1),i-1)*K_Z_U(gen_sUp(flav2),j-1))-
			   K_uI*K_uI*K_Z_R(1,k)*K_Z_R(1,l)*
			   (K_Z_U(gen_sUp(flav1),i-1)*K_Z_U(gen_sUp(flav2),j-1) + 
			    K_Z_U(gen_sUp(flav1)+3,i-1)*K_Z_U(gen_sUp(flav2)+3,j-1)));
	      
	      kcpl1 = kcpl0;
	      
	      vertex[vanz].cpl[0]  = kcpl0; 
	      vertex[vanz].cpl[1]  = kcpl1;
	      vertex[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();
	      
	      
	      vertex[vanz].Color.push_back(Color_Function(cf::D));     
	      vertex[vanz].Color.back().SetParticleArg(1,2);     
	      vertex[vanz].Color.back().SetStringArg('1','2');     
	      
	      
	      vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("SSSS",LF_Key()));     
	      
	      vertex[vanz].on      = 1;
		vertex.push_back(Single_Vertex());vanz++;
	    }
	  }
	}
      }
    }
  }
  //down-squarks Higgs interactions
  for (short int i=1;i<7;++i) {
    if (i<4) s_qu=1000000 + 2*i - 1;
    else     s_qu=2000000 + 2*i - 7;
    Flavour flav1 = Flavour((kf_code)(s_qu));
    for (short int j=i;j<7;++j) {
      if (j<4) s_qu=1000000 + 2*j - 1;
      else     s_qu=2000000 + 2*j - 7;
      Flavour flav2 = Flavour((kf_code)(s_qu));
      if (flav1.IsOn() && flav2.IsOn() && gen_sDown(flav1)==gen_sDown(flav2)) {
		
	//Hmin -> sdown - sdown - Hmin 
	if (flHmin.IsOn()) {
	  vertex[vanz].in[0] = flHmin;
	  vertex[vanz].in[1] = flav1.Bar();
	  vertex[vanz].in[2] = flav2;
	  vertex[vanz].in[3] = flHmin;
	  
	  vertex[vanz].nleg  = 4;  
	  
	  help = K_zero;
	  if(i==j) help = num_1;
	 
	  Kabbala fac1,fac2,fac3;
	  fac1=fac2=fac3=K_zero;
	  
	  for (int x=0;x<3;x++) {
	    for (int y=0;y<3;y++) {
	      fac1 += K_u(0)*conj_K_CKM(0,x)*K_CKM(0,y)*K_Z_D(y,i-1)*K_Z_D(x,j-1); 
	      fac2 += K_u(1)*conj_K_CKM(1,x)*K_CKM(1,y)*K_Z_D(y,i-1)*K_Z_D(x,j-1); 
	      fac3 += K_u(2)*conj_K_CKM(2,x)*K_CKM(2,y)*K_Z_D(y,i-1)*K_Z_D(x,j-1); 
	    }
	  }
	  
	  kcpl0 = M_I*(g1*g1/(costW*costW*num_6)*K_A_H(0,0)*
		       (help-(num_3-num_2*sintW*sintW)/(num_2*sintW*sintW)*
			K_Z_D(gen_sDown(flav1),i-1)*K_Z_D(gen_sDown(flav1),j-1))
		       -K_d(gen_sDown(flav1))*K_d(gen_sDown(flav1))*K_Z_H(1,0)*K_Z_H(1,0)*
		       K_Z_D(gen_sDown(flav1)+3,i-1)*K_Z_D(gen_sDown(flav1)+3,j-1)
		       -K_Z_H(1,0)*K_Z_H(1,0)*
		       (fac1 + fac2 + fac3));
	         
	  kcpl1 = kcpl0;
	  
	  vertex[vanz].cpl[0]  = kcpl0; 
	  vertex[vanz].cpl[1]  = kcpl1;
	  vertex[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();
	  
	  
	  vertex[vanz].Color.push_back(Color_Function(cf::D));     
	  vertex[vanz].Color.back().SetParticleArg(1,2);     
	  vertex[vanz].Color.back().SetStringArg('1','2');     
	  
	  
	  vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("SSSS",LF_Key()));     
	  
	  vertex[vanz].on      = 1;
	  vertex.push_back(Single_Vertex());vanz++;
	}
	//A0 -> sdown - sdown - A0 
	if (flA0.IsOn()) {
	  vertex[vanz].in[0] = flA0;
	  vertex[vanz].in[1] = flav1.Bar();
	  vertex[vanz].in[2] = flav2;
	  vertex[vanz].in[3] = flA0;
	  
	  vertex[vanz].nleg  = 4;  
	  
	  help = K_zero;
	  if(i==j) help = num_1;
	  
	  
	  kcpl0 = M_I*(g1*g1/(costW*costW*num_6)*K_A_H(0,0)*
		       (help+(num_3-num_4*sintW*sintW)/(num_2*sintW*sintW)*
			K_Z_D(gen_sDown(flav1),i-1)*K_Z_D(gen_sDown(flav1),j-1))
		       -K_d(gen_sDown(flav1))*K_d(gen_sDown(flav1))*K_Z_H(1,0)*K_Z_H(1,0)*
		       (K_Z_D(gen_sDown(flav1),i-1)*K_Z_D(gen_sDown(flav1),j-1) + 
			K_Z_D(gen_sDown(flav1)+3,i-1)*K_Z_D(gen_sDown(flav1)+3,j-1)));
	  
	  kcpl1 = kcpl0;
	  
	  vertex[vanz].cpl[0]  = kcpl0; 
	  vertex[vanz].cpl[1]  = kcpl1;
	  vertex[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();
	  
	  
	  vertex[vanz].Color.push_back(Color_Function(cf::D));     
	  vertex[vanz].Color.back().SetParticleArg(1,2);     
	  vertex[vanz].Color.back().SetStringArg('1','2');     
	    
	  
	  vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("SSSS",LF_Key()));     
	  
	  vertex[vanz].on      = 1;
	  vertex.push_back(Single_Vertex());vanz++;
	}
	//h0/H0 -> sdown - sdown - h0/H0 
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
	      
	      kcpl0 = M_I*(g1*g1/(costW*costW*num_6)*K_A_R(k,l)*
			   (help+(num_3-num_4*sintW*sintW)/(num_2*sintW*sintW)*
			    K_Z_D(gen_sDown(flav1),i-1)*K_Z_D(gen_sDown(flav2),j-1))-
			   K_d(gen_sDown(flav1))*K_d(gen_sDown(flav1))*K_Z_R(1,k)*K_Z_R(1,l)*
			   (K_Z_D(gen_sDown(flav1),i-1)*K_Z_D(gen_sDown(flav2),j-1) + 
			    K_Z_D(gen_sDown(flav1)+3,i-1)*K_Z_D(gen_sDown(flav2)+3,j-1)));
	      
	      kcpl1 = kcpl0;
	      
	      vertex[vanz].cpl[0]  = kcpl0; 
	      vertex[vanz].cpl[1]  = kcpl1;
	      vertex[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();
	      
	      
	      vertex[vanz].Color.push_back(Color_Function(cf::D));     
	      vertex[vanz].Color.back().SetParticleArg(1,2);     
	      vertex[vanz].Color.back().SetStringArg('1','2');     
	      
	      
	      vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("SSSS",LF_Key()));     
	      
	      vertex[vanz].on      = 1;
		vertex.push_back(Single_Vertex());vanz++;
	    }
	  }
	}
      }
    }
  }
  //H- -> sdown - supb - H0/h0
  if (flHmin.IsOn()) {
    for (int i=1;i<7;i++) {
      if (i<4) s_qu=1000000 + 2*i;
      else     s_qu=2000000 + 2*i - 6;
      Flavour flsup = Flavour((kf_code)(s_qu));
      for (int j=1;j<7;j++) {
	if (j<4) s_qu=1000000 + 2*j - 1;
	else     s_qu=2000000 + 2*j - 7;
	Flavour flsdown = Flavour((kf_code)(s_qu));
	if (flsup.IsOn() && flsdown.IsOn()) {
	  
	  Kabbala fac1,fac2,fac3,fac4;
	  fac1=fac2=fac3=fac4=K_zero;
	  
	  for (int x=0;x<3;x++) {
	    for (int y=0;y<3;y++) {
	      fac1 += conj_K_CKM(x,y)*K_Z_U(x,i-1)*K_Z_D(y,j-1); 
	      fac2 += conj_K_CKM(x,y)*K_u(x)*K_d(y)*K_Z_U(x+3,i-1)*K_Z_D(y+3,j-1); 
	      fac3 += conj_K_CKM(x,y)*K_u(x)*K_u(x)*K_Z_U(x,i-1)*K_Z_D(y,j-1);
	      fac4 += conj_K_CKM(x,y)*K_d(y)*K_d(y)*K_Z_U(x,i-1)*K_Z_D(y,j-1);
	    }
	  }
	  
	  for (int l=0;l<2;l++) {
	    Flavour flh = Flavour((kf_code)(25+l*10));
	    if (flh.IsOn()) {
	      
	      vertex[vanz].in[0] = flHmin;
	      vertex[vanz].in[1] = flsup.Bar();
	      vertex[vanz].in[2] = flsdown;
	      vertex[vanz].in[3] = flh;
	      
	      vertex[vanz].nleg  = 4;  
	      
	      
	      kcpl0 = M_I/root2*(-g2*g2/num_2*(K_Z_H(0,0)*K_Z_R(0,l) + 
					       K_Z_H(1,0)*K_Z_R(1,l))*fac1 
				 - K_A_P(l,0)*fac2 
				 + K_Z_H(1,0)*K_Z_R(1,l)*fac3
				 + K_Z_H(0,0)*K_Z_R(0,l)*fac4);
	      
	      kcpl1 = kcpl0;
	    
	      vertex[vanz].cpl[0]  = kcpl0; 
	      vertex[vanz].cpl[1]  = kcpl1;
	      vertex[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();
	      
	      
	      vertex[vanz].Color.push_back(Color_Function(cf::D));     
	      vertex[vanz].Color.back().SetParticleArg(1,2);     
	      vertex[vanz].Color.back().SetStringArg('1','2');     
	      
	      
	      vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("SSSS",LF_Key()));     
	      
	      vertex[vanz].on      = 1;
	      vertex.push_back(Single_Vertex());vanz++;
	    }
	  }
	  if (flA0.IsOn()) {
	     vertex[vanz].in[0] = flHmin;
	     vertex[vanz].in[1] = flsup.Bar();
	     vertex[vanz].in[2] = flsdown;
	     vertex[vanz].in[3] = flA0;
	     
	     vertex[vanz].nleg  = 4;  
	     
	      
	     kcpl0 = -num_1/root2*(g2*g2/num_2*K_A_H(0,0)*fac1 
				   + K_Z_H(1,0)*K_Z_R(1,0)*fac3
				   - K_Z_H(0,0)*K_Z_R(0,0)*fac4);
	     
	     kcpl1 = kcpl0;
	     
	     vertex[vanz].cpl[0]  = kcpl0; 
	     vertex[vanz].cpl[1]  = kcpl1;
	     vertex[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();
	     
	     
	     vertex[vanz].Color.push_back(Color_Function(cf::D));     
	     vertex[vanz].Color.back().SetParticleArg(1,2);     
	     vertex[vanz].Color.back().SetStringArg('1','2');     
	     
	     
	     vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("SSSS",LF_Key()));     
	     
	     vertex[vanz].on      = 1;
	     vertex.push_back(Single_Vertex());vanz++;
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
		 ComplexMatrixElement(std::string("CKM"),i,j));
} 
  
Kabbala Interaction_Model_sQuark_EW::conj_K_CKM(short int i,short int j)       
{   
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("V^\\ti_{")+string(hi)+string(hj)+string("}"),
		 conj(ComplexMatrixElement(std::string("CKM"),i,j)));
} 
 

Kabbala Interaction_Model_sQuark_EW::K_Z_D(short int i,short int j)       
{   
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("Z^{")+string(hi)+string(hj)+string("}_D"),
		 ComplexMatrixElement(std::string("Z_d"),i,j));
}  

Kabbala Interaction_Model_sQuark_EW::K_Z_U(short int i,short int j)       
{   
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("Z^{")+string(hi)+string(hj)+string("}_U"),
		 ComplexMatrixElement(std::string("Z_u"),i,j));
}  

Kabbala Interaction_Model_sQuark_EW::K_A_H(short int i,short int j) {
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("A^{")+string(hi)+string(hj)+string("}_H"),
		 ComplexMatrixElement(std::string("Z_H"),0,i) *
		 ComplexMatrixElement(std::string("Z_H"),0,j) -
		 ComplexMatrixElement(std::string("Z_H"),1,i) *
		 ComplexMatrixElement(std::string("Z_H"),1,j));
}  

Kabbala Interaction_Model_sQuark_EW::K_A_R(short int i,short int j) {
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("A^{")+string(hi)+string(hj)+string("}_R"),
		 ComplexMatrixElement(std::string("Z_R"),0,i) *
		 ComplexMatrixElement(std::string("Z_R"),0,j) -
		 ComplexMatrixElement(std::string("Z_R"),1,i) *
		 ComplexMatrixElement(std::string("Z_R"),1,j));
}  

Kabbala Interaction_Model_sQuark_EW::K_A_P(short int i,short int j) {
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("A^{")+string(hi)+string(hj)+string("}_P"),
		 ComplexMatrixElement(std::string("Z_R"),0,i) *
		 ComplexMatrixElement(std::string("Z_H"),1,j) +
		 ComplexMatrixElement(std::string("Z_R"),1,i) *
		 ComplexMatrixElement(std::string("Z_H"),0,j));
}  

Kabbala Interaction_Model_sQuark_EW::K_w_S(short int i,short int j)
{   
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("w^{")+string(hi)+string(hj)+string("}_S"),
		 ComplexMatrixElement(std::string("ws"),i,j));

}  

Kabbala Interaction_Model_sQuark_EW::K_u_S(short int i,short int j)
{   
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("u^{")+string(hi)+string(hj)+string("}_S"),
		 ComplexMatrixElement(std::string("us"),i,j));

}  

Kabbala Interaction_Model_sQuark_EW::K_e_S(short int i,short int j)
{   
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("e^{")+string(hi)+string(hj)+string("}_S"),
		 ComplexMatrixElement(std::string("es"),i,j));
}  

Kabbala Interaction_Model_sQuark_EW::K_d_S(short int i,short int j)
{   
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("d^{")+string(hi)+string(hj)+string("}_S"),
		 ComplexMatrixElement(std::string("ds"),i,j));
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
		 ComplexMatrixElement(std::string("Z_H"),i,j));
}     

Kabbala Interaction_Model_sQuark_EW::K_Z_R(short int i,short int j) {   
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("Z^{")+string(hi)+string(hj)+string("}_R"),
		 ComplexMatrixElement(std::string("Z_R"),i,j));
}  

Kabbala Interaction_Model_sQuark_EW::K_B_R(short int i) {   
  char hi[2];
  sprintf(hi,"%i",i);
  return Kabbala(string("B^{")+string(hi)+string("}_R"),
		 v1.Value() * ComplexMatrixElement(std::string("Z_R"),0,i) -
		 v2.Value() * ComplexMatrixElement(std::string("Z_R"),1,i) );
}  

Kabbala Interaction_Model_sQuark_EW::K_Z_PL(short int i,short int j)       
{   
  char hi[2];
  char hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("Z^\\p_{")+string(hi)+string(hj)+string("}"),
		 ComplexMatrixElement(std::string("Z^+"),i,j));
}  

Kabbala Interaction_Model_sQuark_EW::K_Z_MI(short int i,short int j)       
{   
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("Z^\\m_{")+string(hi)+string(hj)+string("}"),
		 ComplexMatrixElement(std::string("Z^-"),i,j));
}  

Kabbala Interaction_Model_sQuark_EW::K_Z_N(short int i,short int j)       
{   
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("Z^{")+string(hi)+string(hj)+string("}_N"),
		 ComplexMatrixElement(std::string("Z^N"),i,j));
}  
//we use transposed convention !!! 

int Interaction_Model_sQuark_EW::gen_sUp(Flavour fl)
{
  int gen_sUp;
  
  if (fl.Kfcode() == 1000002 || fl.Kfcode() == 2000002)
    gen_sUp = 0;
  if (fl.Kfcode() == 1000004 || fl.Kfcode() == 2000004)
    gen_sUp = 1;
  if (fl.Kfcode() == 1000006 || fl.Kfcode() == 2000006)
    gen_sUp = 2;

  return gen_sUp;
}

int Interaction_Model_sQuark_EW::gen_sDown(Flavour fl)
{
  int gen_sDown;

  if (fl.Kfcode() == 1000001 || fl.Kfcode() == 2000001)
    gen_sDown = 0;
  if (fl.Kfcode() == 1000003 || fl.Kfcode() == 2000003)
    gen_sDown = 1;
  if (fl.Kfcode() == 1000005 || fl.Kfcode() == 2000005)
    gen_sDown = 2;

  return gen_sDown;
}

Kabbala Interaction_Model_sQuark_EW::K_u(short int i)
{
  char hi[2];
  sprintf(hi,"%i",i);
  
  return Kabbala(string("u^")+string(hi),
		 Flavour((kf_code)(2*i+2)).Yuk()/v2.Value()*sqrt(2.));
}

Kabbala Interaction_Model_sQuark_EW::K_d(short int i)
{
  char hi[2];
  sprintf(hi,"%i",i);
  
  return Kabbala(string("d^")+string(hi),
		 -Flavour((kf_code)(2*i+1)).Yuk()/v1.Value()*sqrt(2.));
}
	
