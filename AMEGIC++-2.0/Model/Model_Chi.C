#include "Model_Chi.H"
#include <stdio.h>

using namespace AMEGIC;
using namespace APHYTOOLS;
using namespace std;

// RK + SS Vertices complete

Model_Chi::~Model_Chi()
{
  if (ini_higgs==1) delete moHiggs;
}

void Model_Chi::Init()
{
  if (moHiggs==NULL) {
    //THDM part
    moHiggs = new Model_Higgs;
    moHiggs->Init();
    ini_higgs = 1;
  }
  
  // Spectrum

  if (isa!=NULL) SpCh.Interface(isa);
  else
  SpCh.Init();
  if (isa!=NULL) SpNe.Interface(isa);
  else
  SpNe.Init();

  g1     = Kabbala(string("e"),sqrt(4.*M_PI*moHiggs->Aqed()));
  g2     = Kabbala(string("e/sin\\theta_W"), g1.Value()/K_sinTW().Value());
  PL     = Kabbala(string("P_L"),1.);
  PR     = Kabbala(string("P_R"),1.);
  M_I    = Kabbala(string("i"),Complex(0.,1.));
  root2  = Kabbala(string("\\sqrt{2}"),sqrt(2.));
  K_zero = Kabbala(string("zero"),0.);
  num_2  = Kabbala(string("2"),2.);    	
  num_4  = Kabbala(string("4"),4.);    		
}

void Model_Chi::Init(Model_Higgs* _moHiggs,Isajet* _isa)
{
  moHiggs = _moHiggs;
  isa = _isa;

  Init();
}

void Model_Chi::c_FFS(Single_Vertex* v,int& vanz)
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
	  v[vanz].in[0] = flav2;
	  v[vanz].in[1] = flav3;
	  v[vanz].in[2] = flav1;
	
	  kcpl0 = -M_I/root2*g2*
	    (K_Z_R(0,k-31)*K_Z_PL(0,i-41)*K_Z_MI(1,j-41)+
	     K_Z_R(1,k-31)*K_Z_PL(1,i-41)*K_Z_MI(0,j-41));
	    
	  kcpl1 = -M_I/root2*g2*
	    (K_Z_R(0,k-31)*K_Z_PL(0,j-41)*K_Z_MI(1,i-41)+
	     K_Z_R(1,k-31)*K_Z_PL(1,j-41)*K_Z_MI(0,i-41));
	  
	  v[vanz].cpl[0] = kcpl0.Value();
	  v[vanz].cpl[1] = kcpl1.Value();
	  v[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();
	  v[vanz].cpl[2] = 0.;v[vanz].cpl[3]  = 0.;

	  v[vanz].ncf   = 1;
	  v[vanz].Color = new Color_Function; 
	  
	  v[vanz].Color->type = cf::None; 
	  
	  v[vanz].nlf     = 1;
	  v[vanz].Lorentz = new Lorentz_Function; 
	  
	  v[vanz].Lorentz->type       = lf::FFS;     
	  
	  v[vanz].on     = 1;
	  vanz++;
	  //checked RK
	}
	if (flav1.IsOn() && flav2.IsOn() && flav3.IsOn() && (k == 33)) {
	  v[vanz].in[0] = flav2;
	  v[vanz].in[1] = flav3;
	  v[vanz].in[2] = flav1;

	  kcpl0 = -g2/root2*
	    (K_Z_H(0,0)*K_Z_PL(0,i-41)*K_Z_MI(1,j-41)+
	     K_Z_H(1,0)*K_Z_PL(1,i-41)*K_Z_MI(0,j-41));
	    
	  kcpl1 = -g2/root2*
	    (K_Z_H(0,0)*K_Z_PL(0,j-41)*K_Z_MI(1,i-41)+
	     K_Z_H(1,0)*K_Z_PL(1,j-41)*K_Z_MI(0,i-41));
	  
	  v[vanz].cpl[0] = kcpl0.Value();
	  v[vanz].cpl[1] = kcpl1.Value();
	  v[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();
	  v[vanz].cpl[2] = 0.;v[vanz].cpl[3]  = 0.;
	  
	  v[vanz].ncf   = 1;
	  v[vanz].Color = new Color_Function; 
	  
	  v[vanz].Color->type = cf::None; 
	  
	  v[vanz].nlf     = 1;
	  v[vanz].Lorentz = new Lorentz_Function; 
	  
	  v[vanz].Lorentz->type       = lf::FFS;     
	  
	  v[vanz].on     = 1;
	  vanz++;
	}
      }
    }
    
  // Neutralino - Neutralino - Higgs
  for (short int i=43;i<47;i++) {
    Flavour flav1 = Flavour(kf::code(i));
      for (short int j=43;j<47;j++) {
	Flavour flav2 = Flavour(kf::code(j));
	if (flav1.IsOn() && flav2.IsOn() && flav3.IsOn() && (k != 33) && (i<=j)) {
	  v[vanz].in[0] = flav1;
	  v[vanz].in[2] = flav2;
	  v[vanz].in[1] = flav3;
	  // Z_N have to be complex conjugated for cpl[0]

	  kcpl0 = K_yuk_sign(flav1)*M_I/(K_cosTW()*num_2)*g2*
	    ((K_Z_R(0,k-31)*K_Z_N(2,i-43)-K_Z_R(1,k-31)*K_Z_N(3,i-43))*
	     (K_Z_N(0,j-43)*K_sinTW()-K_Z_N(1,j-43)*K_cosTW())+
	     (K_Z_R(0,k-31)*K_Z_N(2,j-43)-K_Z_R(1,k-31)*K_Z_N(3,j-43))*
	     (K_Z_N(0,i-43)*K_sinTW()-K_Z_N(1,i-43)*K_cosTW()));
	  
	  kcpl1 = K_yuk_sign(flav2)*M_I/(K_cosTW()*num_2)*g2*
	    ((K_Z_R(0,k-31)*K_Z_N(2,i-43)-K_Z_R(1,k-31)*K_Z_N(3,i-43))*
	     (K_Z_N(0,j-43)*K_sinTW()-K_Z_N(1,j-43)*K_cosTW())+
	     (K_Z_R(0,k-31)*K_Z_N(2,j-43)-K_Z_R(1,k-31)*K_Z_N(3,j-43))*
	     (K_Z_N(0,i-43)*K_sinTW()-K_Z_N(1,i-43)*K_cosTW()));
	  		  	  
	  v[vanz].cpl[0] = kcpl0.Value();
	  v[vanz].cpl[1] = kcpl1.Value();
	  v[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();
	  v[vanz].cpl[2] = 0.;v[vanz].cpl[3]  = 0.;
	 
	  v[vanz].ncf   = 1;
	  v[vanz].Color = new Color_Function; 
	  
	  v[vanz].Color->type = cf::None; 
	  
	  v[vanz].nlf     = 1;
	  v[vanz].Lorentz = new Lorentz_Function; 
	  
	  v[vanz].Lorentz->type       = lf::FFS;     
	  
	  v[vanz].on     = 1;
	  vanz++;
	  //checked FK & RK & SS (new -kcpl0)
	}
	if (flav1.IsOn() && flav2.IsOn() && flav3.IsOn() && (k == 33)) {
	  v[vanz].in[0] = flav1;
	  v[vanz].in[2] = flav2;
	  v[vanz].in[1] = flav3;
	  // Z_N have to be complex conjugated for cpl[0]

	  kcpl0 = -K_yuk_sign(flav1)/(K_cosTW()*num_2)*g2*
	    ((K_Z_H(0,k-33)*K_Z_N(2,i-43)-K_Z_H(1,k-33)*K_Z_N(3,i-43))*
	     (K_Z_N(0,j-43)*K_sinTW()-K_Z_N(1,j-43)*K_cosTW()) +
	     (K_Z_H(0,k-33)*K_Z_N(2,j-43)-K_Z_H(1,k-33)*K_Z_N(3,j-43))*
	     (K_Z_N(0,i-43)*K_sinTW()-K_Z_N(1,i-43)*K_cosTW()));

	  kcpl1 = K_yuk_sign(flav2)/(K_cosTW()*num_2)*g2*
	    ((K_Z_H(0,k-33)*K_Z_N(2,i-43)-K_Z_H(1,k-33)*K_Z_N(3,i-43))*
	     (K_Z_N(0,j-43)*K_sinTW()-K_Z_N(1,j-43)*K_cosTW()) +
	     (K_Z_H(0,k-33)*K_Z_N(2,j-43)-K_Z_H(1,k-33)*K_Z_N(3,j-43))*
	     (K_Z_N(0,i-43)*K_sinTW()-K_Z_N(1,i-43)*K_cosTW()));

	  v[vanz].cpl[0] = kcpl0.Value(); 
	  v[vanz].cpl[1] = kcpl1.Value();
	  v[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();
	  v[vanz].cpl[2]  = 0.;v[vanz].cpl[3]  = 0.;
	  
	  v[vanz].ncf   = 1;
	  v[vanz].Color = new Color_Function; 
	  
	  v[vanz].Color->type = cf::None; 
	  
	  v[vanz].nlf     = 1;
	  v[vanz].Lorentz = new Lorentz_Function; 
	  
	  v[vanz].Lorentz->type       = lf::FFS;     
	  
	  v[vanz].on      = 1;
	  vanz++;
	  //checked SS (new -kclp0)
	}
      }
    }
  }
  // Chargino - Higgs - Neutralino
  Flavour flHm = Flavour(kf::Hmin);
  if (flHm.IsOn()) {
    for (short int i=41;i<43;i++) {
      Flavour flav1 = Flavour(kf::code(i));
      for (short int j=43;j<47;j++) {
	Flavour flav2 = Flavour(kf::code(j));
	if (flav1.IsOn() && flav2.IsOn()) {
	  v[vanz].in[0] = flav1;
	  v[vanz].in[1] = flHm;
	  v[vanz].in[2] = flav2;
	  	    
	  kcpl0 = -M_I*g2/K_cosTW()*
	    K_Z_H(1,0)*(K_Z_PL(1,i-41)/root2*
			(K_Z_N(0,j-43)*K_sinTW()+K_Z_N(1,j-43)*K_cosTW())+
			K_Z_PL(0,i-41)*K_Z_N(3,j-43)*K_cosTW());
	 	  
	  kcpl1 = M_I*g2/K_cosTW()*
	    K_Z_H(0,0)*(K_Z_MI(1,i-41)/root2*
			(K_Z_N(0,j-43)*K_sinTW()+K_Z_N(1,j-43)*K_cosTW())-
			K_Z_MI(0,i-41)*K_Z_N(2,j-43)*K_cosTW());
	  
	  v[vanz].cpl[0] = kcpl0.Value(); 
	  v[vanz].cpl[1] = kcpl1.Value();
	  v[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();
	  v[vanz].cpl[2] = 0.;v[vanz].cpl[3]  = 0.;
	  
	  v[vanz].ncf   = 1;
	  v[vanz].Color = new Color_Function; 
	  
	  v[vanz].Color->type = cf::None; 
	  
	  v[vanz].nlf     = 1;
	  v[vanz].Lorentz = new Lorentz_Function; 
	  
	  v[vanz].Lorentz->type       = lf::FFS;     
	  
	  v[vanz].on     = 1;
	  vanz++;
	}
      }
    }
  }
}

void Model_Chi::c_FFV(Single_Vertex* v,int& vanz)
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
	  v[vanz].in[0] = flav1;
	  v[vanz].in[1] = flph;
	  v[vanz].in[2] = flav2;
	  
	  Kabbala charge1 = Kabbala(string("Q_{")+flav1.TexName()+string("}"),
				    flav1.Charge());
	  
	  //fixed with Chargino pair production !!!
	  kcpl0 = -M_I*g1*charge1;
	  kcpl1 = kcpl0;
	  	  
	  v[vanz].cpl[0]  = kcpl0.Value();
	  v[vanz].cpl[1]  = kcpl1.Value();
	  v[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();
	  v[vanz].cpl[2]  = 0.;v[vanz].cpl[3]  = 0.;
	  
	  v[vanz].ncf   = 1;
	  v[vanz].Color = new Color_Function; 
	  
	  v[vanz].Color->type = cf::None; 
	  
	  v[vanz].nlf     = 1;
	  v[vanz].Lorentz = new Lorentz_Function; 
		
	  v[vanz].Lorentz->type       = lf::Gamma;     
	  v[vanz].Lorentz->SetParticleArg(1);     
	  
	  v[vanz].on      = 1;
	  vanz++;
	  //checked RK + SS (tested by chargino pair production)
	}
	if (flZ.IsOn()) {	
	  Kabbala charge1 = Kabbala(string("Q_{")+flav1.TexName()+string("}"),
				  flav1.Charge());
	  Kabbala helper = Kabbala(string("zero"),0.);
	  
	  if (flav1 == flav2) helper = Kabbala(string("cos2\\theta_W"),
					       1.-2.*K_sinTW().Value()*
					       K_sinTW().Value());
	  v[vanz].in[0] = flav1;
	  v[vanz].in[1] = flZ;
	  v[vanz].in[2] = flav2;	  

	  //fixed with Chargino pair production !!!
	    
	  kcpl0 = -M_I/(K_cosTW()*num_2)*
	    g2*charge1*(K_Z_PL(0,j-41)*K_Z_PL(0,i-41) + helper);
	  
	  kcpl1 = -M_I/(K_cosTW()*num_2)*g2*charge1*
	    (K_Z_MI(0,j-41)*K_Z_MI(0,i-41) + helper);
	  	  
	  v[vanz].cpl[0] = kcpl0.Value();
	  v[vanz].cpl[1] = kcpl1.Value(); 
	  v[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();
	  v[vanz].cpl[2] = 0.;v[vanz].cpl[3]  = 0.;
	  
	  v[vanz].ncf   = 1;
	  v[vanz].Color = new Color_Function; 
	  
	  v[vanz].Color->type = cf::None; 
	  
	  v[vanz].nlf     = 1;
	  v[vanz].Lorentz = new Lorentz_Function; 
	  
	  v[vanz].Lorentz->type       = lf::Gamma;     
	  v[vanz].Lorentz->SetParticleArg(1);     
	  
	  v[vanz].on     = 1;
	  vanz++;
	  //checked RK + SS 
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
	v[vanz].in[0] = flav1;
	v[vanz].in[1] = Flavour(kf::W);
	v[vanz].in[2] = flav2.Bar();

	Kabbala charge2 = Kabbala(string("Q_{")+flav2.TexName()+string("}"),
				  flav2.Charge());
	
	kcpl0 = M_I*g2*(K_Z_N(1,i-43)*K_Z_MI(0,j-41)+
				K_Z_N(2,i-43)*K_Z_MI(1,j-41)/root2);
	
	kcpl1 = M_I*g2*(K_Z_N(1,i-43)*K_Z_PL(0,j-41)-
				K_Z_N(3,i-43)*K_Z_PL(1,j-41)/root2);
	
	v[vanz].cpl[0] = kcpl0.Value();
	v[vanz].cpl[1] = kcpl1.Value(); 
	v[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();
	v[vanz].cpl[2] = 0.;v[vanz].cpl[3]  = 0.;
	
	v[vanz].ncf   = 1;
	v[vanz].Color = new Color_Function; 
	  
	v[vanz].Color->type = cf::None; 
	
	v[vanz].nlf     = 1;
	v[vanz].Lorentz = new Lorentz_Function; 
	
	v[vanz].Lorentz->type       = lf::Gamma;     
	v[vanz].Lorentz->SetParticleArg(1);     
	
	v[vanz].on     = 1;
	vanz++;
      }

    }
  }

  //Neutralino - Z - Neutralino
  for (short int j=43;j<47;j++) {
    Flavour flav1 = Flavour(kf::code(j));
    for (short int i=j;i<47;i++) {
      Flavour flav2 = Flavour(kf::code(i));
      if (flav1.IsOn() && flav2.IsOn()) {
		
	v[vanz].in[0] = flav1;
	v[vanz].in[1] = Flavour(kf::Z);	
	v[vanz].in[2] = flav2;
		
	kcpl0 = -M_I/(K_cosTW()*num_2)*g2*(K_Z_N_com(3,i-43)*K_Z_N_com_conj(3,j-43)-
					   K_Z_N_com(2,i-43)*K_Z_N_com_conj(2,j-43)
					   );
	kcpl1 = M_I/(K_cosTW()*num_2)*g2*(K_Z_N_com_conj(3,i-43)*K_Z_N_com(3,j-43)-
					  K_Z_N_com_conj(2,i-43)*K_Z_N_com(2,j-43)
					  );
	v[vanz].cpl[0] = kcpl0.Value();
	v[vanz].cpl[1] = kcpl1.Value();
	v[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();
	v[vanz].cpl[2]  = 0.;v[vanz].cpl[3]  = 0.;
	
	v[vanz].ncf   = 1;
	v[vanz].Color = new Color_Function; 
	  
	v[vanz].Color->type = cf::None; 
	
	v[vanz].nlf     = 1;
	v[vanz].Lorentz = new Lorentz_Function; 
	
	v[vanz].Lorentz->type       = lf::Gamma;     
	v[vanz].Lorentz->SetParticleArg(1);     
	
	v[vanz].on      = 1;
	vanz++;
	//checked FK & RK & SS 
      }
    }
  }  
}


Kabbala Model_Chi::K_Z_PL(short int i,short int j)       
{   
  char hi[2];
  char hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("Z^\\p_{")+string(hi)+string(hj)+string("}"),SpCh.Zplus(i,j));
}  
Kabbala Model_Chi::K_Z_MI(short int i,short int j)       
{   
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("Z^\\m_{")+string(hi)+string(hj)+string("}"),SpCh.Zminus(i,j));
}  
Kabbala Model_Chi::K_Z_N(short int i,short int j)       
{   
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("Z^{")+string(hi)+string(hj)+string("}_N"),SpNe.Z_N(i,j));
}  
//we use transposed convention !!! 

Kabbala Model_Chi::K_Z_N_com(short int i,short int j)       
{   
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);

  //Complex exp_i = (Flavour(kf::code(43+i)).MassSign()==1) ? 1 : Complex(0.,1.);
  //Complex exp_i(cos(M_PI*(1.-etai)/4.),sin(M_PI*(1.-etai)/4.));
  Complex exp_i = Complex(1.,0.);
  
  return Kabbala(string("Z^{")+string(hi)+string(hj)+string("}_N"),exp_i*SpNe.Z_N(i,j));
}  
Kabbala Model_Chi::K_Z_N_com_conj(short int i,short int j)       
{   
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  //Complex exp_i = (Flavour(kf::code(43+i)).MassSign()==1) ? 1 : Complex(0.,-1.);
  //Complex exp_i(cos(M_PI*(1.-etai)/4.),-sin(M_PI*(1.-etai)/4.));
  Complex exp_i = Complex(1.,0.);
  
  return Kabbala(string("Z^{")+string(hi)+string(hj)+string("}_N"),exp_i*SpNe.Z_N(i,j));
}  










