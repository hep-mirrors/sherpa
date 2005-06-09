#include "Interaction_Model_Inos.H"
#include "MathTools.H"
#include "Message.H"
#include "Run_Parameter.H"
#include <stdio.h>

using namespace AMEGIC;
using namespace ATOOLS;
using namespace std;

Interaction_Model_Inos::Interaction_Model_Inos(MODEL::Model_Base * _model,
					       std::string _cplscheme,std::string _yukscheme) :
  Interaction_Model_Base(_model,_cplscheme,_yukscheme)
{ 
  double Ecms2 = sqr(rpa.gen.Ecms());

  g1     = Kabbala(string("g_1"),
		   sqrt(4.*M_PI*ScalarFunction(std::string("alpha_QED"),Ecms2)));
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

  for (short int k=31;k<34;k++) {
    Flavour flav3 = Flavour(kf::code(k));
    // Chargino - Chargino - Higgs
    for (short int i=41;i<43;i++) {
      Flavour flav1 = Flavour(kf::code(i));
      for (short int j=i;j<43;j++) {
	Flavour flav2 = Flavour(kf::code(j));
	if (flav1.IsOn() && flav2.IsOn() && flav3.IsOn() && (k != 33)) {
	  vertex[vanz].in[0] = flav2;
	  vertex[vanz].in[1] = flav3;
	  vertex[vanz].in[2] = flav1;
	
	  kcpl1 = -M_I/root2*g2*
	    (K_Z_R(0,k-31)*K_Z_PL(0,i-41)*K_Z_MI(1,j-41)+
	     K_Z_R(1,k-31)*K_Z_PL(1,i-41)*K_Z_MI(0,j-41));
	    
	  kcpl0 = -M_I/root2*g2*
	    (K_Z_R(0,k-31)*K_Z_PL(0,j-41)*K_Z_MI(1,i-41)+
	     K_Z_R(1,k-31)*K_Z_PL(1,j-41)*K_Z_MI(0,i-41));
	  
	  vertex[vanz].cpl[0] = kcpl0.Value();
	  vertex[vanz].cpl[1] = kcpl1.Value();
	  vertex[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();
	  vertex[vanz].cpl[2] = 0.;vertex[vanz].cpl[3]  = 0.;

	  vertex[vanz].ncf   = 1;
	  vertex[vanz].Color = new Color_Function(cf::None); 
	  
	  vertex[vanz].nlf     = 1;
	  vertex[vanz].Lorentz = new Lorentz_Function(lf::FFS);
	  
	  vertex[vanz].on     = 1;
	  vertex.push_back(Single_Vertex());vanz++;
	  //checked RK
	}
	if (flav1.IsOn() && flav2.IsOn() && flav3.IsOn() && (k == 33)) {
	  vertex[vanz].in[0] = flav2;
	  vertex[vanz].in[1] = flav3;
	  vertex[vanz].in[2] = flav1;

	  kcpl0 = g2/root2*
	    (K_Z_H(0,0)*K_Z_PL(0,i-41)*K_Z_MI(1,j-41)+
	     K_Z_H(1,0)*K_Z_PL(1,i-41)*K_Z_MI(0,j-41));
	    
	  kcpl1 = -g2/root2*
	    (K_Z_H(0,0)*K_Z_PL(0,j-41)*K_Z_MI(1,i-41)+
	     K_Z_H(1,0)*K_Z_PL(1,j-41)*K_Z_MI(0,i-41));
	  
	  vertex[vanz].cpl[0] = kcpl0.Value();
	  vertex[vanz].cpl[1] = kcpl1.Value();
	  vertex[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();
	  vertex[vanz].cpl[2] = 0.;vertex[vanz].cpl[3]  = 0.;
	  
	  vertex[vanz].ncf   = 1;
	  vertex[vanz].Color = new Color_Function(cf::None); 
	  
	  vertex[vanz].nlf     = 1;
	  vertex[vanz].Lorentz = new Lorentz_Function(lf::FFS);
	  
	  vertex[vanz].on     = 1;
	  vertex.push_back(Single_Vertex());vanz++;
	}
      }
    }
    
  // Neutralino - Neutralino - Higgs
  for (short int i=43;i<47;i++) {
    Flavour flav1 = Flavour(kf::code(i));
      for (short int j=43;j<47;j++) {
	Flavour flav2 = Flavour(kf::code(j));
	if (flav1.IsOn() && flav2.IsOn() && flav3.IsOn() && (k != 33) && (i<=j)) {
	  vertex[vanz].in[0] = flav1;
	  vertex[vanz].in[2] = flav2;
	  vertex[vanz].in[1] = flav3;
	  // Z_N have to be complex conjugated for cpl[0]
	  
	  kcpl0 = M_I/(costW*num_2)*g2*
	    ((K_Z_R(0,k-31)*K_Z_N(2,j-43)-K_Z_R(1,k-31)*K_Z_N(3,j-43))*
	     (K_Z_N(0,i-43)*sintW-K_Z_N(1,i-43)*costW)+
	     (K_Z_R(0,k-31)*K_Z_N(2,i-43)-K_Z_R(1,k-31)*K_Z_N(3,i-43))*
	     (K_Z_N(0,j-43)*sintW-K_Z_N(1,j-43)*costW));
	  
	  kcpl1 = M_I/(costW*num_2)*g2*
	    ((K_Z_R(0,k-31)*K_Z_N(2,i-43)-K_Z_R(1,k-31)*K_Z_N(3,i-43))*
	     (K_Z_N(0,j-43)*sintW-K_Z_N(1,j-43)*costW)+
	     (K_Z_R(0,k-31)*K_Z_N(2,j-43)-K_Z_R(1,k-31)*K_Z_N(3,j-43))*
	     (K_Z_N(0,i-43)*sintW-K_Z_N(1,i-43)*costW));
	  		  	  
	  vertex[vanz].cpl[0] = kcpl0.Value();
	  vertex[vanz].cpl[1] = kcpl1.Value();
	  vertex[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();
	  vertex[vanz].cpl[2] = 0.;vertex[vanz].cpl[3]  = 0.;
	 
	  vertex[vanz].ncf   = 1;
	  vertex[vanz].Color = new Color_Function(cf::None); 
	  
	  vertex[vanz].nlf     = 1;
	  vertex[vanz].Lorentz = new Lorentz_Function(lf::FFS);
	  
	  vertex[vanz].on     = 1;
	  vertex.push_back(Single_Vertex());vanz++;
	  //checked FK & RK & SS 
	}
	if (flav1.IsOn() && flav2.IsOn() && flav3.IsOn() && (k == 33)) {
	  vertex[vanz].in[0] = flav1;
	  vertex[vanz].in[2] = flav2;
	  vertex[vanz].in[1] = flav3;
	  // Z_N have to be complex conjugated for cpl[0]

	  kcpl0 = -g2/(costW*num_2)*
	    ((K_Z_H(0,k-33)*K_Z_N(2,j-43)-K_Z_H(1,k-33)*K_Z_N(3,j-43))*
	     (K_Z_N(0,i-43)*sintW-K_Z_N(1,i-43)*costW) +
	     (K_Z_H(0,k-33)*K_Z_N(2,i-43)-K_Z_H(1,k-33)*K_Z_N(3,i-43))*
	     (K_Z_N(0,j-43)*sintW-K_Z_N(1,j-43)*costW));

	  kcpl1 = g2/(costW*num_2)*
	    ((K_Z_H(0,k-33)*K_Z_N(2,i-43)-K_Z_H(1,k-33)*K_Z_N(3,i-43))*
	     (K_Z_N(0,j-43)*sintW-K_Z_N(1,j-43)*costW) +
	     (K_Z_H(0,k-33)*K_Z_N(2,j-43)-K_Z_H(1,k-33)*K_Z_N(3,j-43))*
	     (K_Z_N(0,i-43)*sintW-K_Z_N(1,i-43)*costW));

	  vertex[vanz].cpl[0] = kcpl0.Value(); 
	  vertex[vanz].cpl[1] = kcpl1.Value();
	  vertex[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();
	  vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;
	  
	  vertex[vanz].ncf   = 1;
	  vertex[vanz].Color = new Color_Function(cf::None); 
	  
	  vertex[vanz].nlf     = 1;
	  vertex[vanz].Lorentz = new Lorentz_Function(lf::FFS);
	  
	  vertex[vanz].on      = 1;
	  vertex.push_back(Single_Vertex());vanz++;
	  //checked SS 
	}
      }
    }
  }
  // Neutralino - Higgs - Chargino
  Flavour flHm = Flavour(kf::Hmin);
  if (flHm.IsOn()) {
    for (short int i=41;i<43;i++) {
      Flavour flav1 = Flavour(kf::code(i));
      for (short int j=43;j<47;j++) {
	Flavour flav2 = Flavour(kf::code(j));
	if (flav1.IsOn() && flav2.IsOn()) {
	  vertex[vanz].in[0] = flav2;
	  vertex[vanz].in[1] = flHm.Bar();
	  vertex[vanz].in[2] = flav1;
	  
	  kcpl0 = -M_I*g2/costW*
	    K_Z_H(1,0)*(K_Z_PL(1,i-41)/root2*
			(K_Z_N(0,j-43)*sintW+K_Z_N(1,j-43)*costW)+
			K_Z_PL(0,i-41)*K_Z_N(3,j-43)*costW);
	  
	  kcpl1 = M_I*g2/costW*
	    K_Z_H(0,0)*(K_Z_MI(1,i-41)/root2*
			(K_Z_N(0,j-43)*sintW+K_Z_N(1,j-43)*costW)-
			K_Z_MI(0,i-41)*K_Z_N(2,j-43)*costW);
	  
	  vertex[vanz].cpl[0] = kcpl0.Value(); 
	  vertex[vanz].cpl[1] = kcpl1.Value();
	  vertex[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();
	  vertex[vanz].cpl[2] = 0.;vertex[vanz].cpl[3]  = 0.;
	  
	  vertex[vanz].ncf   = 1;
	  vertex[vanz].Color = new Color_Function(cf::None); 
	  
	  vertex[vanz].nlf     = 1;
	  vertex[vanz].Lorentz = new Lorentz_Function(lf::FFS);
	  
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
   Flavour flph = Flavour(kf::photon);
   Flavour flZ = Flavour(kf::Z);

   //Chargino - Z/Photon - Chargino
   for (short int i=41;i<43;i++) {
     Flavour flav1 = Flavour(kf::code(i));
     for (short int j=i;j<43;j++) {
       Flavour flav2 = Flavour(kf::code(j));
       if (flav1.IsOn() && flav2.IsOn()) {
	 if (flav1==flav2 && flph.IsOn()) {
	   vertex[vanz].in[0] = flav1;
	   vertex[vanz].in[1] = flph;
	   vertex[vanz].in[2] = flav2;
	   
	   kcpl0 = M_I*g1;
	   kcpl1 = kcpl0;
	   
	   vertex[vanz].cpl[0]  = kcpl0.Value();
	   vertex[vanz].cpl[1]  = kcpl1.Value();
	   vertex[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();
	   vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;
	   
	   vertex[vanz].ncf   = 1;
	   vertex[vanz].Color = new Color_Function(cf::None); 
	   
	   vertex[vanz].nlf     = 1;
	   vertex[vanz].Lorentz = new Lorentz_Function(lf::Gamma);
	   vertex[vanz].Lorentz->SetParticleArg(1);     
	   
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
	   
	   kcpl1 = M_I/(costW*num_2)*g2*
	     (K_Z_MI(0,j-41)*K_Z_MI(0,i-41) + helper);

	   kcpl0 = M_I/(costW*num_2)*g2*
	     (K_Z_PL(0,j-41)*K_Z_PL(0,i-41) + helper);
	   
	   vertex[vanz].cpl[0] = kcpl0.Value();
	   vertex[vanz].cpl[1] = kcpl1.Value(); 
	   vertex[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();
	   vertex[vanz].cpl[2] = 0.;vertex[vanz].cpl[3]  = 0.;
	   
	   vertex[vanz].ncf   = 1;
	   vertex[vanz].Color = new Color_Function(cf::None); 
	   
	   vertex[vanz].nlf     = 1;
	   vertex[vanz].Lorentz = new Lorentz_Function(lf::Gamma);
	   vertex[vanz].Lorentz->SetParticleArg(1);     
	   
	   vertex[vanz].on     = 1;
	   vertex.push_back(Single_Vertex());vanz++;
	 }     
       } 
     }
   }
   
   //Chargino - Neutralino - W-
   for (short int i=43;i<47;i++) {
     Flavour flav1 = Flavour(kf::code(i));
     for (short int j=41;j<43;j++) {
       Flavour flav2 = Flavour(kf::code(j));
       if (flav1.IsOn() && flav2.IsOn()) {
	 vertex[vanz].in[0] = flav1;
	 vertex[vanz].in[1] = Flavour(kf::W).Bar();
	 vertex[vanz].in[2] = flav2;

	 kcpl0 = M_I*g2*(K_Z_N(1,i-43)*K_Z_MI(0,j-41)+
			 K_Z_N(2,i-43)*K_Z_MI(1,j-41)/root2);
	 
	 kcpl1 = M_I*g2*(K_Z_N(1,i-43)*K_Z_PL(0,j-41)-
			 K_Z_N(3,i-43)*K_Z_PL(1,j-41)/root2);

	 vertex[vanz].cpl[0] = kcpl0.Value();
	 vertex[vanz].cpl[1] = kcpl1.Value(); 
	 vertex[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();
	 vertex[vanz].cpl[2] = 0.;vertex[vanz].cpl[3]  = 0.;
	 
	 vertex[vanz].ncf   = 1;
	 vertex[vanz].Color = new Color_Function(cf::None); 
	 
	 vertex[vanz].nlf     = 1;
	 vertex[vanz].Lorentz = new Lorentz_Function(lf::Gamma);
	 vertex[vanz].Lorentz->SetParticleArg(1);     
	 
	 vertex[vanz].on     = 1;
	 vertex.push_back(Single_Vertex());vanz++;
       }
       
     }
   }
   
   //Neutralino - Z - Neutralino
   for (short int j=43;j<47;j++) {
     Flavour flav1 = Flavour(kf::code(j));
     for (short int i=j;i<47;i++) {
       Flavour flav2 = Flavour(kf::code(i));
       if (flav1.IsOn() && flav2.IsOn() && flZ.IsOn()) {
	 
	 vertex[vanz].in[0] = flav1;
	 vertex[vanz].in[1] = Flavour(kf::Z);	
	 vertex[vanz].in[2] = flav2;
	 
	 kcpl0 = -M_I/(costW*num_2)*g2*
	   (K_Z_N_com(3,i-43)*K_Z_N_com_conj(3,j-43)-
	    K_Z_N_com(2,i-43)*K_Z_N_com_conj(2,j-43));
	 kcpl1 = M_I/(costW*num_2)*g2*
	   (K_Z_N_com_conj(3,i-43)*K_Z_N_com(3,j-43)-
	    K_Z_N_com_conj(2,i-43)*K_Z_N_com(2,j-43));
	 
	 vertex[vanz].cpl[0] = kcpl0.Value();
	 vertex[vanz].cpl[1] = kcpl1.Value();
	 vertex[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();
	 vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;
	 
	 vertex[vanz].ncf   = 1;
	 vertex[vanz].Color = new Color_Function(cf::None); 
	 
	 vertex[vanz].nlf     = 1;
	 vertex[vanz].Lorentz = new Lorentz_Function(lf::Gamma);
	 vertex[vanz].Lorentz->SetParticleArg(1);     
	 
	 vertex[vanz].on      = 1;
	 vertex.push_back(Single_Vertex());vanz++;
	 //checked FK & RK & SS 
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






