#include "MODEL/Interaction_Models/Interaction_Model_Inos.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include <stdio.h>

using namespace MODEL;
using namespace ATOOLS;
using namespace std;

Interaction_Model_Inos::Interaction_Model_Inos(MODEL::Model_Base * _model,
					       std::string _cplscheme,std::string _yukscheme) :
  Interaction_Model_Base("",_model,_cplscheme,_yukscheme)
{ 
  g1     = Kabbala(string("g_1"),
		   sqrt(4.*M_PI*ScalarFunction(std::string("alpha_QED"),rpa.gen.CplScale())));
  g2     = Kabbala(string("g_1/\\sin\\theta_W"), 
		   g1.Value()/sqrt(ScalarConstant(std::string("sin2_thetaW"))));
  sintW  = Kabbala(std::string("\\sin\\theta_W"),
		   sqrt(ScalarConstant(std::string("sin2_thetaW"))));
  costW  = Kabbala(std::string("\\cos\\theta_W"),
		   sqrt(1.-ScalarConstant(std::string("sin2_thetaW"))));
  PL     = Kabbala(string("P_L"),1.);
  PR     = Kabbala(string("P_R"),1.);
  M_I    = Kabbala(string("i"),Complex(0.,1.));
  root2  = Kabbala(string("\\sqrt{2}"),sqrt(2.));
  K_zero = Kabbala(string("zero"),0.);
  num_2  = Kabbala(string("2"),2.);    	
  num_4  = Kabbala(string("4"),4.);    		
}

void Interaction_Model_Inos::c_FFS(std::vector<Single_Vertex>& vertex,int& vanz)
{
  Kabbala kcpl0,kcpl1;
  int higgs,n_ino,c_ino;

  for (short int k=0;k<3;k++) {
    if (k==0) higgs=25;
    if (k==1) higgs=35;
    if (k==2) higgs=36;
    Flavour flav3 = Flavour((kf_code)(higgs));
    // Chargino - Chargino - Higgs
    for (short int i=0;i<2;i++) {
      if (i==0) c_ino=1000024;
      else      c_ino=1000037;
      Flavour flav1 = Flavour((kf_code)(c_ino));
      for (short int j=i;j<2;j++) {
	if (j==0) c_ino=1000024;
	else      c_ino=1000037;
      	Flavour flav2 = Flavour((kf_code)(c_ino));
	if (flav1.IsOn() && flav2.IsOn() && flav3.IsOn() && (k!=2)) {
	  vertex[vanz].in[0] = flav1;
	  vertex[vanz].in[1] = flav3;
	  vertex[vanz].in[2] = flav2;
	
	  kcpl0 = -M_I/root2*g2*
	    (K_Z_R(0,k)*K_Z_PL(0,j)*K_Z_MI(1,i)+
	     K_Z_R(1,k)*K_Z_PL(1,j)*K_Z_MI(0,i));
	  
	  kcpl1 = -M_I/root2*g2*
	    (K_Z_R(0,k)*K_Z_PL(0,i)*K_Z_MI(1,j)+
	     K_Z_R(1,k)*K_Z_PL(1,i)*K_Z_MI(0,j));
	
	  vertex[vanz].cpl[0] = kcpl0;
	  vertex[vanz].cpl[1] = kcpl1;
	  vertex[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();

	  
	  vertex[vanz].Color.push_back(Color_Function(cf::None)); 
	  
	  
	  vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("FFS",LF_Key()));
	  
	  vertex[vanz].on     = 1;
	  vertex.push_back(Single_Vertex());vanz++;
	}
	if (flav1.IsOn() && flav2.IsOn() && flav3.IsOn() && (k==2)) {
	  vertex[vanz].in[0] = flav1;
	  vertex[vanz].in[2] = flav3;
	  vertex[vanz].in[1] = flav2;

	  kcpl0 = g2/root2*
	    (K_Z_H(0,0)*K_Z_PL(0,j)*K_Z_MI(1,i)+
	     K_Z_H(1,0)*K_Z_PL(1,j)*K_Z_MI(0,i));

	  kcpl1 = -g2/root2*
	    (K_Z_H(0,0)*K_Z_PL(0,i)*K_Z_MI(1,j)+
	     K_Z_H(1,0)*K_Z_PL(1,i)*K_Z_MI(0,j));
	  
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
    
  // Neutralino - Neutralino - Higgs
  for (short int i=0;i<4;i++) {
    if (i<2)  n_ino=1000022+i;
    if (i==2) n_ino=1000025;
    if (i==3) n_ino=1000035;
    Flavour flav1 = Flavour((kf_code)(n_ino));
      for (short int j=i;j<4;j++) {
	if (j<2)  n_ino=1000022+j;
	if (j==2) n_ino=1000025;
	if (j==3) n_ino=1000035;
    	Flavour flav2 = Flavour((kf_code)(n_ino));
	if (flav1.IsOn() && flav2.IsOn() && flav3.IsOn() && (k!=2)) {
	  vertex[vanz].in[0] = flav1;
	  vertex[vanz].in[2] = flav2;
	  vertex[vanz].in[1] = flav3;
	  // Z_N have to be complex conjugated for cpl[0]
	  
	  kcpl0 = M_I/(costW*num_2)*g2*
	    ((K_Z_R(0,k)*K_Z_N(2,j)-K_Z_R(1,k)*K_Z_N(3,j))*
	     (K_Z_N(0,i)*sintW-K_Z_N(1,i)*costW)+
	     (K_Z_R(0,k)*K_Z_N(2,i)-K_Z_R(1,k)*K_Z_N(3,i))*
	     (K_Z_N(0,j)*sintW-K_Z_N(1,j)*costW));
	  
	  kcpl1 = M_I/(costW*num_2)*g2*
	    ((K_Z_R(0,k)*K_Z_N(2,i)-K_Z_R(1,k)*K_Z_N(3,i))*
	     (K_Z_N(0,j)*sintW-K_Z_N(1,j)*costW)+
	     (K_Z_R(0,k)*K_Z_N(2,j)-K_Z_R(1,k)*K_Z_N(3,j))*
	     (K_Z_N(0,i)*sintW-K_Z_N(1,i)*costW));

	  vertex[vanz].cpl[0] = kcpl0;
	  vertex[vanz].cpl[1] = kcpl1;
	  vertex[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();
	 
	  
	  vertex[vanz].Color.push_back(Color_Function(cf::None)); 
	  
	  
	  vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("FFS",LF_Key()));
	  
	  vertex[vanz].on     = 1;
	  vertex.push_back(Single_Vertex());vanz++;
	}
	if (flav1.IsOn() && flav2.IsOn() && flav3.IsOn() && (k==2)) {
	  vertex[vanz].in[0] = flav1;
	  vertex[vanz].in[2] = flav2;
	  vertex[vanz].in[1] = flav3;
	  // Z_N have to be complex conjugated for cpl[0]

	  kcpl0 = -g2/(costW*num_2)*
	    ((K_Z_H(0,0)*K_Z_N(2,j)-K_Z_H(1,0)*K_Z_N(3,j))*
	     (K_Z_N(0,i)*sintW-K_Z_N(1,i)*costW) +
	     (K_Z_H(0,0)*K_Z_N(2,i)-K_Z_H(1,0)*K_Z_N(3,i))*
	     (K_Z_N(0,j)*sintW-K_Z_N(1,j)*costW));

	  kcpl1 = g2/(costW*num_2)*
	    ((K_Z_H(0,0)*K_Z_N(2,i)-K_Z_H(1,0)*K_Z_N(3,i))*
	     (K_Z_N(0,j)*sintW-K_Z_N(1,j)*costW) +
	     (K_Z_H(0,0)*K_Z_N(2,j)-K_Z_H(1,0)*K_Z_N(3,j))*
	     (K_Z_N(0,i)*sintW-K_Z_N(1,i)*costW));

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
  // Neutralino - Higgs - Chargino
  Flavour flHm = Flavour(kf_Hplus).Bar();
  if (flHm.IsOn()) {
    for (short int i=0;i<2;i++) {
      if (i==0) c_ino=1000024;
      else      c_ino=1000037;
      Flavour flav1 = Flavour((kf_code)(c_ino));
      for (short int j=0;j<4;j++) {
	if (j<2)  n_ino=1000022+j;
	if (j==2) n_ino=1000025;
	if (j==3) n_ino=1000035;
	Flavour flav2 = Flavour((kf_code)(n_ino));
	if (flav1.IsOn() && flav2.IsOn()) {
	  vertex[vanz].in[0] = flav2;
	  vertex[vanz].in[1] = flHm;
	  vertex[vanz].in[2] = flav1;
	  
	  kcpl0 = -M_I*g2/costW*
	    K_Z_H(1,0)*(K_Z_PL(1,i)/root2*
			(K_Z_N(0,j)*sintW+K_Z_N(1,j)*costW)+
			K_Z_PL(0,i)*K_Z_N(3,j)*costW);
	  
	  kcpl1 = M_I*g2/costW*
	    K_Z_H(0,0)*(K_Z_MI(1,i)/root2*
			(K_Z_N(0,j)*sintW+K_Z_N(1,j)*costW)-
			K_Z_MI(0,i)*K_Z_N(2,j)*costW);

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

void Interaction_Model_Inos::c_FFV(std::vector<Single_Vertex>& vertex,int& vanz)
{
   Kabbala kcpl0,kcpl1;  
   Flavour flph = Flavour(kf_photon);
   Flavour flZ = Flavour(kf_Z);
   int c_ino,n_ino;
   
   //Chargino - Z/Photon - Chargino
   for (short int i=0;i<2;i++) {
     if (i==0) c_ino = 1000024;
     else      c_ino = 1000037;
     Flavour flav1 = Flavour((kf_code)(c_ino));
     for (short int j=i;j<2;j++) {
     if (j==0) c_ino = 1000024;
     else      c_ino = 1000037;
     Flavour flav2 = Flavour((kf_code)(c_ino));
       if (flav1.IsOn() && flav2.IsOn()) {
	 if (flav1==flav2 && flph.IsOn()) {
	   vertex[vanz].in[0] = flav1;
	   vertex[vanz].in[1] = flph;
	   vertex[vanz].in[2] = flav2;
	   
	   kcpl0 = -M_I*g1;
	   kcpl1 = kcpl0;
	   
	   vertex[vanz].cpl[0]  = kcpl0;
	   vertex[vanz].cpl[1]  = kcpl1;
	   vertex[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();
	   
	   
	   vertex[vanz].Color.push_back(Color_Function(cf::None)); 
	   
	   
	   vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("Gamma",LF_Key()));
	   vertex[vanz].Lorentz.back()->SetParticleArg(1);     
	   
	   vertex[vanz].on      = 1;
	   vertex.push_back(Single_Vertex());vanz++;
	 }
	 if (flZ.IsOn()) {	

	   Kabbala helper = Kabbala(string("zero"),0.);
	   
	   if (flav1 == flav2) helper = Kabbala(string("cos2\\theta_W"),
						1.-2.*sintW.Value()*
						sintW.Value());
	   
	   vertex[vanz].in[0] = flav1;
	   vertex[vanz].in[1] = flZ;
	   vertex[vanz].in[2] = flav2;	  

	   kcpl0 = -M_I/(costW*num_2)*g2*
	     (K_Z_MI(0,j)*K_Z_MI(0,i) + helper);

	   kcpl1 = -M_I/(costW*num_2)*g2*
	     (K_Z_PL(0,j)*K_Z_PL(0,i) + helper);
	
	   vertex[vanz].cpl[0] = kcpl0;
	   vertex[vanz].cpl[1] = kcpl1; 
	   vertex[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();
	   
	   
	   vertex[vanz].Color.push_back(Color_Function(cf::None)); 
	   
	   
	   vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("Gamma",LF_Key()));
	   vertex[vanz].Lorentz.back()->SetParticleArg(1);     
	   
	   vertex[vanz].on     = 1;
	   vertex.push_back(Single_Vertex());vanz++;
	 }     
       } 
     }
   }
   
   //Chargino - Neutralino - W+
   for (short int i=0;i<4;i++) {
     if (i<2)  n_ino=1000022+i;
     if (i==2) n_ino=1000025;
     if (i==3) n_ino=1000035;
     Flavour flav1 = Flavour((kf_code)(n_ino));
     for (short int j=0;j<2;j++) {
       if (j==0) c_ino = 1000024;
       else      c_ino = 1000037;
       Flavour flav2 = Flavour((kf_code)(c_ino));
       if (flav1.IsOn() && flav2.IsOn()) {
	 vertex[vanz].in[0] = flav1;
	 vertex[vanz].in[1] = Flavour(kf_Wplus).Bar();
	 vertex[vanz].in[2] = flav2;

	 kcpl0 = M_I*g2*(K_Z_N(1,i)*K_Z_MI(0,j)+
			 K_Z_N(2,i)*K_Z_MI(1,j)/root2);

	 kcpl1 = M_I*g2*(K_Z_N(1,i)*K_Z_PL(0,j)-
			 K_Z_N(3,i)*K_Z_PL(1,j)/root2);
	 
	 vertex[vanz].cpl[0] = kcpl0;
	 vertex[vanz].cpl[1] = kcpl1; 
	 vertex[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();
	 
	 
	 vertex[vanz].Color.push_back(Color_Function(cf::None)); 
	 
	 
	 vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("Gamma",LF_Key()));
	 vertex[vanz].Lorentz.back()->SetParticleArg(1);     
	 
	 vertex[vanz].on     = 1;
	 vertex.push_back(Single_Vertex());vanz++;
       }
       
     }
   }
   
   //Neutralino - Z - Neutralino
   for (short int j=0;j<4;j++) {
     if (j<2)  n_ino=1000022+j;
     if (j==2) n_ino=1000025;
     if (j==3) n_ino=1000035;
     Flavour flav1 = Flavour((kf_code)(n_ino));
     for (short int i=j;i<4;i++) {
       if (i<2)  n_ino=1000022+i;
       if (i==2) n_ino=1000025;
       if (i==3) n_ino=1000035;
       Flavour flav2 = Flavour((kf_code)(n_ino));
       if (flav1.IsOn() && flav2.IsOn() && flZ.IsOn()) {
	 
	 vertex[vanz].in[0] = flav1;
	 vertex[vanz].in[1] = Flavour(kf_Z);	
	 vertex[vanz].in[2] = flav2;
	 
	 kcpl0 = -M_I/(costW*num_2)*g2*
	   (K_Z_N_com(3,j)*K_Z_N_com_conj(3,i)-
	    K_Z_N_com(2,j)*K_Z_N_com_conj(2,i));
	 kcpl1 = M_I/(costW*num_2)*g2*
	   (K_Z_N_com_conj(3,j)*K_Z_N_com(3,i)-
	    K_Z_N_com_conj(2,j)*K_Z_N_com(2,i));
	 
	 vertex[vanz].cpl[0] = kcpl0;
	 vertex[vanz].cpl[1] = kcpl1;
	 vertex[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();
	 
	 
	 vertex[vanz].Color.push_back(Color_Function(cf::None)); 
	 
	 
	 vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("Gamma",LF_Key()));
	 vertex[vanz].Lorentz.back()->SetParticleArg(1);     
	 
	 vertex[vanz].on      = 1;
	 vertex.push_back(Single_Vertex());vanz++;
       }
     }
   }  
}




Kabbala Interaction_Model_Inos::K_CKM(short int i,short int j)       
{   
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("V_{")+string(hi)+string(hj)+string("}"),
		 ComplexMatrixElement(string("CKM"),i,j));
} 
  
Kabbala Interaction_Model_Inos::conj_K_CKM(short int i,short int j)       
{   
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("V^\\ti_{")+string(hi)+string(hj)+string("}"),
		 conj(ComplexMatrixElement(string("CKM"),i,j)));
} 
 

Kabbala Interaction_Model_Inos::K_yuk(Flavour fl) {
  return Kabbala(string("M_{"+fl.TexName()+"}"),fl.Yuk());
}

Kabbala Interaction_Model_Inos::K_yuk_sign(Flavour fl) {
  char hi[3];
  sprintf(hi,"%i",fl.MassSign());
  return Kabbala(string(hi),fl.MassSign());
}

Kabbala Interaction_Model_Inos::K_Z_H(short int i,short int j) {   
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("Z^{")+string(hi)+string(hj)+string("}_H"),
		 ComplexMatrixElement(string("Z_H"),i,j));
}     

Kabbala Interaction_Model_Inos::K_Z_R(short int i,short int j) {   
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("Z^{")+string(hi)+string(hj)+string("}_R"),
		 ComplexMatrixElement(string("Z_R"),i,j));
}  

Kabbala Interaction_Model_Inos::K_Z_PL(short int i,short int j)       
{   
  char hi[2];
  char hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("Z^\\p_{")+string(hi)+string(hj)+string("}"),
		 ComplexMatrixElement(string("Z^+"),i,j));
}  

Kabbala Interaction_Model_Inos::K_Z_MI(short int i,short int j)       
{   
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("Z^\\m_{")+string(hi)+string(hj)+string("}"),
		 ComplexMatrixElement(string("Z^-"),i,j));
}  

Kabbala Interaction_Model_Inos::K_Z_N(short int i,short int j)       
{   
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("Z^{")+string(hi)+string(hj)+string("}_N"),
		 ComplexMatrixElement(string("Z^N"),i,j));
}  
//we use transposed convention !!! 

Kabbala Interaction_Model_Inos::K_Z_N_com(short int i,short int j)       
{   
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);

  Complex exp_i = Complex(1.,0.);
  
  return Kabbala(string("Z^{")+string(hi)+string(hj)+string("}_N"),
		 exp_i*ComplexMatrixElement(string("Z^N"),i,j));
}  

Kabbala Interaction_Model_Inos::K_Z_N_com_conj(short int i,short int j)       
{   
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  Complex exp_i = Complex(1.,0.);
  
  return Kabbala(string("Z^{")+string(hi)+string(hj)+string("}_N"),
		 exp_i*ComplexMatrixElement(string("Z^N"),i,j));
}  






