#include "Model_sLep_Chi.H"
//#include <stdio.h>

using namespace AMEGIC;
using namespace APHYTOOLS;
using namespace std;

//Vertices complete SS

Model_sLep_Chi::~Model_sLep_Chi()
{
  if (ini_sm==1) delete mosm;
  if (ini_ch==1) delete moch;
  if (ini_sL==1) delete mosL;
}

void Model_sLep_Chi::Init()
{
  if (mosm==NULL) {
    //SM part
    mosm = new Model_SM;
    mosm->Init();
    ini_sm = 1;
  }
  if (moch==NULL) {
    //Chi part
    moch = new Model_Chi;
    moch->Init();
    ini_ch = 1;
  }
  if (mosL==NULL) {
    //sLepton part
    mosL = new Model_sLeptons;
    mosL->Init();
    ini_sL = 1;
  }

  g1     = Kabbala(string("e"),sqrt(4.*M_PI*mosm->Aqed()));
  g2     = Kabbala(string("e/sin\\theta_W"), g1.Value()/mosm->SinTW());
  PL     = Kabbala(string("P_L"),1.);
  PR     = Kabbala(string("P_R"),1.);
  M_I    = Kabbala(string("i"),Complex(0.,1.));
  root2  = Kabbala(string("\\sqrt{2}"),sqrt(2.));
  K_zero = Kabbala(string("zero"),0.);
}

void Model_sLep_Chi::Init(Model_SM* _mosm,Model_Chi* _moch,Model_sLeptons* _mosL)
{
  mosm = _mosm;
  moch = _moch;
  mosL = _mosL;

  Init();
}

void Model_sLep_Chi::c_FFS(Single_Vertex* v,int& vanz){

  Kabbala kcpl0,kcpl1;
  
  Kabbala K_l1 = Kabbala(string("l^1"),
			 -Flavour(kf::code(11)).Yuk()/(mosL->K_v1()).Value()*sqrt(2.));

  Kabbala K_l2 = Kabbala(string("l^2"),
			 -Flavour(kf::code(13)).Yuk()/(mosL->K_v1()).Value()*sqrt(2.));
  
  Kabbala K_l3 = Kabbala(string("l^3"),
			 -Flavour(kf::code(15)).Yuk()/(mosL->K_v1()).Value()*sqrt(2.));
  
  //neutrino - sneutrino - neutralino
  
  for (short int i=12;i<17;i+=2) {
    Flavour flav1 = Flavour(kf::code(i));
    for (short int j=43;j<47;j++) {
      Flavour flav2 = Flavour(kf::code(j));
      for (short int k=81;k<84;k++) {
	Flavour flav3 = Flavour(kf::code(k));
	if (flav1.IsOn() && flav2.IsOn() && flav3.IsOn()) {
	  v[vanz].in[0] = flav1;
	  v[vanz].in[1] = flav3;
	  v[vanz].in[2] = flav2;
	  
	  kcpl0 = K_zero;
	  kcpl1 = M_I*g2/(K_cosTW()*root2)*K_Z_nue((i-12)/2,k-81)*
	    (K_Z_N(0,j-43)*K_sinTW()-K_Z_N(1,j-43)*K_cosTW());
	  
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
	  //checked SS
	}
      }
    }
  }  
    
  //neutrino - slepton   - Chargino  

  for (short int i=12;i<17;i+=2) {
    Flavour flav1 = Flavour(kf::code(i));
    
    Kabbala K_lI = Kabbala(string("\\frac{\\m M_{")+flav1.TexName()+string("}}{ v_1}\\sqrt{2}"),
			   -Flavour(kf::code(i-1)).Yuk()/(mosL->K_v1()).Value()*sqrt(2.));
   
    for (short int j=41;j<43;j++) {
      Flavour flav2 = Flavour(kf::code(j));
      for (short int k=71;k<77;k++) {
	Flavour flav3 = Flavour(kf::code(k));
	if (flav1.IsOn() && flav2.IsOn() && flav3.IsOn()) {
	  v[vanz].in[0] = flav1;
	  v[vanz].in[1] = flav3;
	  v[vanz].in[2] = flav2.Bar();
	 
	  kcpl0 = K_zero;
	  kcpl1 = K_yuk_sign(flav2)*(-M_I)*(K_Z_MI(0,j-41)*K_Z_L((i-12)/2,k-71)*g2+
					    K_Z_MI(1,j-41)*K_Z_L((i-12)/2+3,k-71)*K_lI);
	 
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
	  //checked SS
	}
      }
    }
  }

  //lepton - Chargino - sneutrino  
 
  for (short int i=11;i<16;i+=2) {
    Flavour flav1 = Flavour(kf::code(i));
  
    Kabbala K_lI = Kabbala(string("\\frac{(\\m M_{")+flav1.TexName()+string(")}}{ v_1}\\sqrt{2}"),
      -flav1.Yuk()/(mosL->K_v1()).Value()*sqrt(2.));
       
    for (short int j=41;j<43;j++) {
      Flavour flav2 = Flavour(kf::code(j));
      for (short int k=81;k<84;k++) {
	Flavour flav3 = Flavour(kf::code(k));
	if (flav1.IsOn() && flav2.IsOn() && flav3.IsOn()) {
	  v[vanz].in[0] = flav1;
	  v[vanz].in[1] = flav3;
	  v[vanz].in[2] = flav2;

	  kcpl0 = -M_I*K_lI*K_Z_MI(1,j-41)*K_Z_nue((i-11)/2,k-81);
	  kcpl1 = -M_I*g2*K_Z_PL(0,j-41)*K_Z_nue((i-11)/2,k-81); 

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
	  v[vanz].in[0] = flav1;
	  v[vanz].in[1] = flav3;
	  v[vanz].in[2] = flav2;

	  Kabbala K_lI = Kabbala(string("\\frac{\\m M_{")+flav1.TexName()+string("}}{ v_1}\\sqrt{2}"),
				     -flav1.Yuk()/(mosL->K_v1()).Value()*sqrt(2.));
	  
	  kcpl0 = -M_I*(g1/K_cosTW()*root2*K_Z_L((i-11)/2+3,k-71)*K_Z_N_com_conj(0,j-43)
			-K_lI*K_Z_L((i-11)/2,k-71)*K_Z_N_com_conj(2,j-43));
	  
	  kcpl1 = M_I*(g2/(K_cosTW()*root2)*K_Z_L((i-11)/2,k-71)*
		       (K_Z_N_com(0,j-43)*K_sinTW()+K_Z_N_com(1,j-43)*K_cosTW())
		       +K_lI*K_Z_L((i-11)/2+3,k-71)*K_Z_N_com(2,j-43));
	  
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
	  //checked RK & SS
	}
      }
    }
  }
}



















