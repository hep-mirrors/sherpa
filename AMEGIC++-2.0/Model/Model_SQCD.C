#include "Model_SQCD.H"
#include <stdio.h>

using namespace AMEGIC;
using namespace APHYTOOLS;
using namespace std;

Model_SQCD::~Model_SQCD()
{
  if (ini_higgs==1) delete moHiggs;
}

void Model_SQCD::Init()
{
  if (moHiggs==NULL) {
    //THDM part
    moHiggs = new Model_Higgs;
    moHiggs->Init();
    ini_higgs = 1;
  }
  // Spectrum
  if (isa!=NULL) SpsD.Interface(isa);
  else
  SpsD.Init();
  if (isa!=NULL) SpsU.Interface(isa);
  else
  SpsU.Init();

  g1     = Kabbala(string("e"),sqrt(4.*M_PI*moHiggs->Aqed()));
  g2     = Kabbala(string("e/sin\\theta_W"), g1.Value()/moHiggs->SinTW());
  g3     = g3 = Kabbala(string("g_3"),sqrt(4.*M_PI*moHiggs->Aqcd()));
  PL     = Kabbala(string("P_L"),1.);
  PR     = Kabbala(string("P_R"),1.);
  M_I    = Kabbala(string("i"),Complex(0.,1.));
  root2  = Kabbala(string("\\sqrt{2}"),sqrt(2.));
  inv_root2 = Kabbala(string("\\frac{1}{\\sqrt{2}}"),1./sqrt(2.));
  K_zero = Kabbala(string("zero"),0.);
  num_2  = Kabbala(string("2"),2.);    	
  num_3  = Kabbala(string("3"),3.);    		
}

void Model_SQCD::Init(Model_Higgs* _moHiggs,Isajet* _isa)
{
  moHiggs = _moHiggs;
  isa = _isa;
  Init();

}

void Model_SQCD::c_FFS(Single_Vertex* v,int& vanz)
{

  Kabbala kcpl0,kcpl1;

  //quark - squark - gluino

  Flavour flgluino = Flavour(kf::code(47));
  if (flgluino.IsOn()) {
  //uptype - sup - gluino
    for (short int i=2;i<7;i+=2) {
      Flavour flav1 = Flavour(kf::code(i));
      for (short int j=51;j<57;j++) {
	Flavour flav2 = Flavour(kf::code(j));
	if (flav1.IsOn() && flav2.IsOn()) {
	v[vanz].in[0] = flav1;
	v[vanz].in[1] = flav2;
	v[vanz].in[2] = flgluino;
	
	kcpl0 = M_I*g3*root2*K_Z_U((i-2)/2+3,j-51);
	kcpl1 = -M_I*g3*root2*K_Z_U((i-2)/2,j-51);
	
	v[vanz].cpl[0]  = kcpl0.Value();
	v[vanz].cpl[1]  = kcpl1.Value();
	v[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	v[vanz].cpl[2]  = 0.;v[vanz].cpl[3]  = 0.;

	v[vanz].ncf   = 1;
	v[vanz].Color = new Color_Function; 
	
	v[vanz].Color->type       = cf::T;     
	v[vanz].Color->SetParticleArg(2,1,0);     
	v[vanz].Color->SetStringArg('2','1','0');     
	
	v[vanz].nlf     = 1;
	v[vanz].Lorentz = new Lorentz_Function; 
	
	v[vanz].Lorentz->type       = lf::FFS;     
	
	v[vanz].on      = 1;
	vanz++;
      }  
      }
    }
  //downtype - sdown - gluino
    for (short int i=1;i<6;i+=2) {
      Flavour flav1 = Flavour(kf::code(i));
      for (short int j=61;j<67;j++) {
	Flavour flav2 = Flavour(kf::code(j));
	if (flav1.IsOn() && flav2.IsOn()) {
	v[vanz].in[0] = flav1;
	v[vanz].in[1] = flav2;
	v[vanz].in[2] = flgluino;
	
	kcpl0 = M_I*g3*root2*K_Z_D((i-1)/2+3,j-61);
	kcpl1 = -M_I*g3*root2*K_Z_D((i-1)/2,j-61);
	
	v[vanz].ncf   = 1;
	v[vanz].Color = new Color_Function; 
	
	v[vanz].Color->type       = cf::T;     
	v[vanz].Color->SetParticleArg(2,1,0);     
	v[vanz].Color->SetStringArg('2','1','0');     

	v[vanz].cpl[0]  = kcpl0.Value();
	v[vanz].cpl[1]  = kcpl1.Value();
	v[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	v[vanz].cpl[2]  = 0.;v[vanz].cpl[3]  = 0.;
	
	v[vanz].nlf     = 1;
	v[vanz].Lorentz = new Lorentz_Function; 
	
	v[vanz].Lorentz->type       = lf::FFS;     
	
	v[vanz].on      = 1;
	vanz++;
	}  
      }
    }
  }
}

void Model_SQCD::c_FFV(Single_Vertex* v,int& vanz)
{
  Kabbala kcpl0,kcpl1;
  
  //gluino - gluon - gluino
  Flavour flgluino = Flavour(kf::code(47));
  if (flgluino.IsOn()) {      
    Flavour flgluon = Flavour(kf::gluon);
    if (flgluon.IsOn()) {
      v[vanz].in[0] = flgluino;
      v[vanz].in[1] = flgluon;
      v[vanz].in[2] = flgluino;
      
      kcpl0 = -g3; 
      kcpl1 = -g3;
      
      v[vanz].cpl[0]  = kcpl0.Value();
      v[vanz].cpl[1]  = kcpl1.Value();
      v[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
      v[vanz].cpl[2]  = 0.;v[vanz].cpl[3]  = 0.;
      
      v[vanz].ncf   = 1;
      v[vanz].Color = new Color_Function; 
      
      v[vanz].Color->type       = cf::F;     
      v[vanz].Color->SetParticleArg(0,1,2);     
      v[vanz].Color->SetStringArg('0','1','2');     
      
      v[vanz].nlf     = 1;
      v[vanz].Lorentz = new Lorentz_Function; 
      
      v[vanz].Lorentz->type       = lf::Gamma;     
      v[vanz].Lorentz->SetParticleArg(1);     
      
      v[vanz].on      = 1;
      vanz++;
    }   
  }
}

void Model_SQCD::c_SSV(Single_Vertex* v,int& vanz)
{
  Kabbala kcpl0,kcpl1;
 
 //squark - Photon - squark

  Flavour flph = Flavour(kf::photon);
  if (flph.IsOn()) {
    //sUpypes
    for (short int i=51 ;i<57;i++) {
      Flavour flav = Flavour(kf::code(i));
      if (flav.IsOn()) {
	v[vanz].in[0] = flav;
	v[vanz].in[1] = flph;
	v[vanz].in[2] = flav;
	
	Kabbala charge = Kabbala(string("Q_{"+flav.TexName()+"}"),flav.Charge());
	
	//changed sign - -> +
	kcpl0 = M_I*charge*g1;
	kcpl1 = kcpl0;
	
	v[vanz].cpl[0]  = kcpl0.Value();
	v[vanz].cpl[1]  = kcpl1.Value();
	v[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	v[vanz].cpl[2]  = 0.;v[vanz].cpl[3]  = 0.;
	
	v[vanz].ncf   = 1;
	v[vanz].Color = new Color_Function; 
	
	v[vanz].Color->type       = cf::D;     
	v[vanz].Color->SetParticleArg(0,2);     
	v[vanz].Color->SetStringArg('0','2');     
	  
	v[vanz].nlf     = 1;
	v[vanz].Lorentz = new Lorentz_Function; 
	
	v[vanz].Lorentz->type       = lf::SSV;     
	v[vanz].Lorentz->SetParticleArg(0,2,1);     
	
	v[vanz].on      = 1;
	vanz++;
      }  
    }
    //sDowntypes
    for (short int i=61 ;i<67;i++) {
      Flavour flav = Flavour(kf::code(i));
      if (flav.IsOn()) {
	v[vanz].in[0] = flav;
	v[vanz].in[1] = flph;
	v[vanz].in[2] = flav;
	
	Kabbala charge = Kabbala(string("Q_{"+flav.TexName()+"}"),flav.Charge());

	kcpl0 = -M_I*charge*g1;
	kcpl1 = kcpl0;
	
	v[vanz].cpl[0]  = kcpl0.Value();
	v[vanz].cpl[1]  = kcpl1.Value();
	v[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	v[vanz].cpl[2]  = 0.;v[vanz].cpl[3]  = 0.;
	
	v[vanz].ncf   = 1;
	v[vanz].Color = new Color_Function; 
	
	v[vanz].Color->type       = cf::D;     
	v[vanz].Color->SetParticleArg(0,2);     
	v[vanz].Color->SetStringArg('0','2');     
	
	v[vanz].nlf     = 1;
	v[vanz].Lorentz = new Lorentz_Function; 
	
	v[vanz].Lorentz->type       = lf::SSV;     
	v[vanz].Lorentz->SetParticleArg(0,2,1);     
	
	v[vanz].on      = 1;
	vanz++;
      }
    }
  }   

  //squark - Z - squark

  Flavour flZ = Flavour(kf::Z);
  if (flZ.IsOn()) {
    //sUpypes
    for (short int i=51;i<57;i++) {
      Flavour flav1 = Flavour(kf::code(i));
      for (short int j=i;j<57;j++) {
	Flavour flav2 = Flavour(kf::code(j));
	if (flav1.IsOn() && flav2.IsOn()) {
	  
	  v[vanz].in[0] = flav1;
	  v[vanz].in[1] = flZ;
	  v[vanz].in[2] = flav2;
	  
	  Kabbala help = K_zero;
	  
	  if (i==j) {help = K_sinTW()*K_sinTW()*num_2/num_3;}  

	  kcpl0 = M_I*g2/K_cosTW()*
	    /*((K_Z_U(0,j-51)*K_Z_U(0,i-51)+
	      K_Z_U(1,j-51)*K_Z_U(1,i-51)+
	      K_Z_U(2,j-51)*K_Z_U(2,i-51))*/ (K_Z_U(gen_sUp(flav2),j-51)*K_Z_U(gen_sUp(flav2),i-51)/num_2-help);
	  kcpl1 = kcpl0;
	  
	  v[vanz].cpl[0]  = kcpl0.Value();
	  v[vanz].cpl[1]  = kcpl1.Value();
	  v[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	  v[vanz].cpl[2]  = 0.;v[vanz].cpl[3]  = 0.;
	
	  v[vanz].ncf   = 1;
	  v[vanz].Color = new Color_Function; 
	  
	  v[vanz].Color->type       = cf::D;     
	  v[vanz].Color->SetParticleArg(0,2);     
	  v[vanz].Color->SetStringArg('0','2');     
	  
	  v[vanz].nlf     = 1;
	  v[vanz].Lorentz = new Lorentz_Function; 
	  
	  v[vanz].Lorentz->type       = lf::SSV;     
	  v[vanz].Lorentz->SetParticleArg(0,2,1);     
	  
	  v[vanz].on      = 1;
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
	  
	  v[vanz].in[0] = flav1;
	  v[vanz].in[1] = flZ;
	  v[vanz].in[2] = flav2;
	  
	  Kabbala help = K_zero;
	  
	  if (i==j) {help = K_sinTW()*K_sinTW()/num_3;}  
	  
	  kcpl0 = M_I*g2/K_cosTW()*
	    /*((K_Z_D(0,j-61)*K_Z_D(0,i-61)+
	      K_Z_D(1,j-61)*K_Z_D(1,i-61)+
	      K_Z_D(2,j-61)*K_Z_D(2,i-61))*/ (K_Z_D(gen_sDown(flav2),j-61)*K_Z_D(gen_sDown(flav2),i-61)/num_2-help);
	  kcpl1 = kcpl0;

	  v[vanz].cpl[0]  = kcpl0.Value();
	  v[vanz].cpl[1]  = kcpl1.Value();
	  v[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	  v[vanz].cpl[2]  = 0.;v[vanz].cpl[3]  = 0.;
	  
	v[vanz].ncf   = 1;
	v[vanz].Color = new Color_Function; 
	
	v[vanz].Color->type       = cf::D;     
	v[vanz].Color->SetParticleArg(0,2);     
	v[vanz].Color->SetStringArg('0','2');     
	
	v[vanz].nlf     = 1;
	v[vanz].Lorentz = new Lorentz_Function; 
	
	v[vanz].Lorentz->type       = lf::SSV;     
	v[vanz].Lorentz->SetParticleArg(0,2,1);     
	
	v[vanz].on      = 1;
	vanz++;
	}  
      }
    }
  }    
  
//check the summing Convention !!!!
  //supquarks - W - sdownquarks
  Flavour flW = Flavour(kf::W);
  if (flW.IsOn()) {
    for (short int i=51;i<57;i++) {
      Flavour flav1 = Flavour(kf::code(i));
      for (short int j=61;j<67;j++) {
	Flavour flav2 = Flavour(kf::code(j));
	if (flav1.IsOn() && flav2.IsOn()) {
	v[vanz].in[0] = flav2;
	v[vanz].in[1] = flW;
	v[vanz].in[2] = flav1;
	
	/*
	  Kabbala help = Kabbala(string("zero"),0.);
	  
	  for (short int t=0;t<3;t++){
	  for (short int z=0;z<3;z++) { 
	  help += K_Z_D(t,j-61)*K_Z_U(z,i-51)*K_CKM(t,z);   
	  }
	}
	//  kcpl0 = -M_I*g2*inv_root2*help;
	*/
	
	kcpl0 = -M_I*g2*inv_root2*(K_Z_D(gen_sDown(flav2),j-61)*K_Z_U(gen_sUp(flav1),i-51))*
	  K_CKM(gen_sDown(flav2),gen_sUp(flav1));
		
	kcpl1 = kcpl0;
	
	v[vanz].cpl[0]  = kcpl0.Value();
	v[vanz].cpl[1]  = kcpl1.Value();
	v[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	v[vanz].cpl[2]  = 0.;v[vanz].cpl[3]  = 0.;
	
	v[vanz].ncf   = 1;
	v[vanz].Color = new Color_Function; 
	
	v[vanz].Color->type       = cf::D;     
	v[vanz].Color->SetParticleArg(0,2);     
	v[vanz].Color->SetStringArg('0','2');     
	  
	v[vanz].nlf     = 1;
	v[vanz].Lorentz = new Lorentz_Function; 
	
	v[vanz].Lorentz->type       = lf::SSV;     
	v[vanz].Lorentz->SetParticleArg(0,2,1);     
	
	v[vanz].on      = 1;
	vanz++;
	} 
    }			
  }
} 

  //sQuark - Gluon - sQuark
   
  Flavour flgl = Flavour(kf::gluon); 
  if (flgl.IsOn()) {

  //uptypes 
  for (short int i=51;i<57;i++) {
    Flavour flav = Flavour(kf::code(i));
    if (flav.IsOn()) { 
      v[vanz].in[0] = flav;
      v[vanz].in[1] = flgl;
      v[vanz].in[2] = flav;
            
      kcpl0 = g3*M_I;
      kcpl1 = kcpl0;
      
      v[vanz].cpl[0]  = kcpl0.Value();
      v[vanz].cpl[1]  = kcpl1.Value();
      v[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
      v[vanz].cpl[2]  = 0.;v[vanz].cpl[3]  = 0.;
      
      v[vanz].ncf   = 1;
      v[vanz].Color = new Color_Function; 
      
      v[vanz].Color->type       = cf::T;     
      v[vanz].Color->SetParticleArg(1,0,2);     
      v[vanz].Color->SetStringArg('1','0','2');     

      v[vanz].nlf     = 1;
      v[vanz].Lorentz = new Lorentz_Function; 
      
      v[vanz].Lorentz->type       = lf::SSV;     
      v[vanz].Lorentz->SetParticleArg(0,2,1);     
      
      v[vanz].on      = 1;
      vanz++;
      
    } 
  }
  //downtypes 
  for (short int i=61;i<67;i++) {
    Flavour flav = Flavour(kf::code(i));
    if (flav.IsOn()) { 
      v[vanz].in[0] = flav;
      v[vanz].in[1] = flgl;
      v[vanz].in[2] = flav;
            
      kcpl0 = -g3*M_I;
      kcpl1 = kcpl0;
      
      v[vanz].cpl[0]  = kcpl0.Value();
      v[vanz].cpl[1]  = kcpl1.Value();
      v[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
      v[vanz].cpl[2]  = 0.;v[vanz].cpl[3]  = 0.;
      
      v[vanz].ncf   = 1;
      v[vanz].Color = new Color_Function; 
      
      v[vanz].Color->type       = cf::T;     
      v[vanz].Color->SetParticleArg(1,0,2);     
      v[vanz].Color->SetStringArg('1','0','2');     

      v[vanz].nlf     = 1;
      v[vanz].Lorentz = new Lorentz_Function; 
      
      v[vanz].Lorentz->type       = lf::SSV;     
      v[vanz].Lorentz->SetParticleArg(0,2,1);     
      
      v[vanz].on      = 1;
      vanz++;
    } 
  }
  }
}

void Model_SQCD::c_SSS(Single_Vertex* v,int& vanz)
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
				 (K_v2()).Value()*sqrt(2.));

	  v[vanz].in[0] = flav1;
	  v[vanz].in[1] = flA0;
	  v[vanz].in[2] = flav2;

	  kcpl0 = -(K_uI*K_Z_H(0,0)*(K_h()*K_Z_U(gen_sUp(flav1),i-51)*K_Z_U(gen_sUp(flav1)+3,j-51)-
				     conj_K_h()*K_Z_U(gen_sUp(flav1),j-51)*
				     K_Z_U(gen_sUp(flav1)+3,i-51))+
		    K_Z_H(1,0)*(K_u_S(gen_sUp(flav1),gen_sUp(flav2))*K_Z_U(gen_sUp(flav1),j-51)*
				K_Z_U(gen_sUp(flav2)+3,i-51)-K_u_S(gen_sUp(flav1),gen_sUp(flav2))*
				K_Z_U(gen_sUp(flav1),i-51)*K_Z_U(gen_sUp(flav2)+3,j-51))+
		    K_Z_H(0,0)*(K_w_S(gen_sUp(flav1),gen_sUp(flav2))*K_Z_U(gen_sUp(flav1),i-51)*
				K_Z_U(gen_sUp(flav2)+3,j-51)-K_w_S(gen_sUp(flav1),gen_sUp(flav2))*
				K_Z_U(gen_sUp(flav1),j-51)*K_Z_U(gen_sUp(flav2)+3,i-51)))
	    *inv_root2;
	  
	  kcpl1 = kcpl0;
	  
	  v[vanz].cpl[0]  = kcpl0.Value();
	  v[vanz].cpl[1]  = kcpl1.Value();
	  v[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	  v[vanz].cpl[2]  = 0.;v[vanz].cpl[3]  = 0.;

	  v[vanz].ncf   = 1;
	  v[vanz].Color = new Color_Function; 
      	
	  v[vanz].Color->type       = cf::D;     
	  v[vanz].Color->SetParticleArg(0,2);     
	  v[vanz].Color->SetStringArg('0','2');     
	  
	  v[vanz].nlf     = 1;
	  v[vanz].Lorentz = new Lorentz_Function; 
	  
	  v[vanz].Lorentz->type       = lf::SSS;     
	
	  v[vanz].on      = 1;
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
				-Flavour(kf::code(2*gen_sDown(flav1)+1)).Yuk()/(K_v1()).Value()*sqrt(2.));
	  
	  v[vanz].in[0] = flav1;
	  v[vanz].in[1] = flA0;
	  v[vanz].in[2] = flav2;

	  kcpl0 = -(K_dI*K_Z_H(1,0)*(conj_K_h()*K_Z_D(gen_sDown(flav1),i-61)*
				     K_Z_D(gen_sDown(flav1)+3,j-61)-
				     K_h()*K_Z_D(gen_sDown(flav1),i-61)*
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
				K_Z_D(gen_sDown(flav2)+3,j-61)))*inv_root2;
	  
	  kcpl1 = kcpl0;

	  v[vanz].cpl[0]  = kcpl0.Value();
	  v[vanz].cpl[1]  = kcpl1.Value();
	  v[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	  v[vanz].cpl[2]  = 0.;v[vanz].cpl[3]  = 0.;
	  
	  v[vanz].ncf   = 1;
	  v[vanz].Color = new Color_Function; 
      	  
	  v[vanz].Color->type       = cf::D;     
	  v[vanz].Color->SetParticleArg(0,2);     
	  v[vanz].Color->SetStringArg('0','2');     
	  
	  v[vanz].nlf     = 1;
	  v[vanz].Lorentz = new Lorentz_Function; 
	  
	  v[vanz].Lorentz->type       = lf::SSS;     
	
	  v[vanz].on      = 1;
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
	  v[vanz].in[0] = flav1;
	  v[vanz].in[1] = flH;
	  v[vanz].in[2] = flav2;

	  Kabbala help = K_zero;

	  Kabbala K_uI = Kabbala(string("u^I"),Flavour(kf::code(2*gen_sUp(flav1)+2)).Yuk()/
				 (K_v2()).Value()*sqrt(2.));

	  Kabbala fac = Kabbala(string("\\frac{3-8sin^2\\theta_W}{4sin^2\\theta_W}"),
				(3.-8.*(K_sinTW()).Value()*(K_sinTW()).Value())/
				(4.*(K_sinTW()).Value()*(K_sinTW()).Value()));

	  if (i==j) {help = Kabbala(string("1"),1.);}

	  kcpl0 = M_I*(-g1*g1/(K_cosTW()*K_cosTW()*num_3)*K_B_R(k-31)*
		       (help+fac*K_Z_U(gen_sUp(flav1),i-51)*K_Z_U(gen_sUp(flav1),j-51))-
		       K_uI*K_uI*K_v2()*K_Z_R(1,k-31)*
		       (K_Z_U(gen_sUp(flav1),i-51)*K_Z_U(gen_sUp(flav1),j-51)+
			K_Z_U(gen_sUp(flav1)+3,i-51)*K_Z_U(gen_sUp(flav1)+3,j-51))+
		       K_Z_R(1,k-31)*inv_root2*K_u_S(gen_sUp(flav1),gen_sUp(flav2))*
						  (K_Z_U(gen_sUp(flav1),i-51)*
						   K_Z_U(gen_sUp(flav2)+3,j-51)+
						   K_Z_U(gen_sUp(flav1),j-51)*
						   K_Z_U(gen_sUp(flav2)+3,i-51))+
		       K_Z_R(0,k-31)*inv_root2*K_w_S(gen_sUp(flav1),gen_sUp(flav2))*
						  (K_Z_U(gen_sUp(flav1),i-51)*
						   K_Z_U(gen_sUp(flav2)+3,j-51)+
						   K_Z_U(gen_sUp(flav1),j-51)*
						   K_Z_U(gen_sUp(flav2)+3,i-51))+
		       K_uI*K_Z_R(0,k-31)*inv_root2*(conj_K_h()*K_Z_U(gen_sUp(flav1),j-51)*
						       K_Z_U(gen_sUp(flav1)+3,i-51)+
						       K_h()*K_Z_U(gen_sUp(flav1),i-51)*
						       K_Z_U(gen_sUp(flav2)+3,j-51)));

	  kcpl1 = kcpl0;
	  	  
	  v[vanz].cpl[0]  = kcpl0.Value();
	  v[vanz].cpl[1]  = kcpl1.Value();
	  v[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();
	  v[vanz].cpl[2]  = 0.;v[vanz].cpl[3]  = 0.;
	
	  v[vanz].ncf   = 1;
	  v[vanz].Color = new Color_Function; 
      	  
	  v[vanz].Color->type       = cf::D;     
	  v[vanz].Color->SetParticleArg(0,2);     
	  v[vanz].Color->SetStringArg('0','2');     
	  
	  v[vanz].nlf     = 1;
	  v[vanz].Lorentz = new Lorentz_Function; 
	  
	  v[vanz].Lorentz->type       = lf::SSS;     
	  
	  v[vanz].on      = 1;
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
	  v[vanz].in[0] = flav1;
	  v[vanz].in[1] = flH;
	  v[vanz].in[2] = flav2;

	  Kabbala help = K_zero;

	  Kabbala K_dI = Kabbala(string("d^I"),
				 -Flavour(kf::code(2*gen_sDown(flav1)+1)).Yuk()/(K_v1()).Value()*sqrt(2.));
	  
	  Kabbala fac = Kabbala(string("\\frac{3-4sin^2\\theta_W}{2sin^2\\theta_W}"),
				(3.-4.*(K_sinTW()).Value()*(K_sinTW()).Value())/
				(2.*(K_sinTW()).Value()*(K_sinTW()).Value()));

	  if (i==j) {help = Kabbala(string("1"),1.);}

	  kcpl0 = M_I*(-g1*g1/(K_cosTW()*K_cosTW()*num_6)*K_B_R(k-31)*
		       (help+fac*K_Z_D(gen_sDown(flav1),i-61)*K_Z_D(gen_sDown(flav1),j-61))-
		       K_dI*K_dI*K_v1()*K_Z_R(0,k-31)*
		       (K_Z_D(gen_sDown(flav1),i-61)*K_Z_D(gen_sDown(flav1),j-61)+
			K_Z_D(gen_sDown(flav1)+3,i-61)*K_Z_D(gen_sDown(flav1)+3,j-61))+
		       K_Z_R(0,k-31)*inv_root2*K_d_S(gen_sDown(flav1),gen_sDown(flav2))*
						  (K_Z_D(gen_sDown(flav1),j-61)*
						   K_Z_D(gen_sDown(flav2)+3,i-61)+
						   K_Z_D(gen_sDown(flav1),i-61)*
						   K_Z_D(gen_sDown(flav2)+3,j-61))+
		       K_Z_R(1,k-31)*inv_root2*K_e_S(gen_sDown(flav1),gen_sDown(flav2))*
						  (K_Z_D(gen_sDown(flav1),j-61)*
						   K_Z_D(gen_sDown(flav2)+3,i-61)+
						   K_Z_D(gen_sDown(flav1),i-61)*
						   K_Z_D(gen_sDown(flav2)+3,j-61))-
		       K_dI*K_Z_R(1,k-31)*inv_root2*(conj_K_h()*K_Z_D(gen_sDown(flav1),i-61)*
						       K_Z_D(gen_sDown(flav1)+3,j-61)+
						       K_h()*K_Z_D(gen_sDown(flav1),j-61)*
						       K_Z_D(gen_sDown(flav2)+3,i-61)));
	  
	  kcpl1 = kcpl0;
	  	  
	  v[vanz].cpl[0]  = kcpl0.Value();
	  v[vanz].cpl[1]  = kcpl1.Value();
	  v[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	  v[vanz].cpl[2]  = 0.;v[vanz].cpl[3]  = 0.;
	  
	  v[vanz].ncf   = 1;
	  v[vanz].Color = new Color_Function; 
      	  
	  v[vanz].Color->type       = cf::D;     
	  v[vanz].Color->SetParticleArg(0,2);     
	  v[vanz].Color->SetStringArg('0','2');     
	  
	  v[vanz].nlf     = 1;
	  v[vanz].Lorentz = new Lorentz_Function; 
	  
	  v[vanz].Lorentz->type       = lf::SSS;     
	  
	  v[vanz].on      = 1;
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
	v[vanz].in[0] = flav1.Bar();
	v[vanz].in[1] = flHmin.Bar();
	v[vanz].in[2] = flav2;
	
	Kabbala K_dI = Kabbala(string("d^I"),
				 -Flavour(kf::code(2*gen_sUp(flav1)+1)).Yuk()/(K_v1()).Value()*sqrt(2.));
		
	Kabbala K_uJ = Kabbala(string("u^I"),Flavour(kf::code(2*gen_sDown(flav2)+2)).Yuk()/
			 (K_v2()).Value()*sqrt(2.));
		
	Kabbala K_massW = Kabbala(string("M_W"),g1.Value()/2.*sqrt(K_v1().Value()*K_v1().Value()+
								   K_v2().Value()*K_v2().Value()));
	
	kcpl0 = M_I*((-(g2*g2)/num_2*(K_v1()*K_Z_H(0,0)+K_v2()*K_Z_H(1,0))+
		      K_v1()*K_dI*K_dI*K_Z_H(0,0)+K_v2()*K_uJ*K_uJ*K_Z_H(1,0))*inv_root2*
		     conj_K_CKM(gen_sUp(flav1),gen_sDown(flav2))*
		     K_Z_D(gen_sUp(flav1),j-61)*K_Z_U(gen_sDown(flav2),i-51)-
		     K_sinTW()*K_massW*root2/g1*K_uJ*K_dI*
		     conj_K_CKM(gen_sUp(flav1),gen_sDown(flav2))*
		     K_Z_D(gen_sUp(flav1)+3,j-61)*K_Z_U(gen_sDown(flav2)+3,i-51)+
		     
		     (K_Z_H(0,0)*conj_K_h()*K_uJ*conj_K_CKM(gen_sUp(flav1),gen_sDown(flav2))+
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
		     K_Z_H(1,0)*K_h()*K_dI*conj_K_CKM(gen_sUp(flav1),gen_sDown(flav2)))*
		     K_Z_U(gen_sDown(flav2),i-51)*K_Z_D(gen_sUp(flav1)+3,j-61));
		     
	  kcpl1 = kcpl0;

	  v[vanz].cpl[0]  = kcpl0.Value();
	  v[vanz].cpl[1]  = kcpl1.Value();
	  v[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	  v[vanz].cpl[2]  = 0.;v[vanz].cpl[3]  = 0.;
	  
	  v[vanz].ncf   = 1;
	  v[vanz].Color = new Color_Function; 
      
	  v[vanz].Color->type       = cf::D;     
	  v[vanz].Color->SetParticleArg(0,2);     
	  v[vanz].Color->SetStringArg('0','2');     
	  
	  v[vanz].nlf     = 1;
	  v[vanz].Lorentz = new Lorentz_Function; 
	  
	  v[vanz].Lorentz->type       = lf::SSS;     
	  
	  v[vanz].on      = 1;
	  vanz++;
      }
    }
  }
}

void Model_SQCD::c_SSVV(Single_Vertex* v,int& vanz)
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
	    
	    v[vanz].in[0] = flavW;
	    v[vanz].in[1] = flav1.Bar();
	    v[vanz].in[2] = flav2;
	    v[vanz].in[3] = flavPhoton;
	    
	    v[vanz].nleg     = 4;
	    
	    kcpl0 = num_2;
	    kcpl1 = kcpl0;
	    
	    v[vanz].cpl[0]  = kcpl0.Value();
	    v[vanz].cpl[1]  = kcpl1.Value();
	    v[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	    v[vanz].cpl[2]  = 0.;v[vanz].cpl[3]  = 0.;
	    
	    v[vanz].ncf   = 1;
	    v[vanz].Color = new Color_Function; 
	    
	    v[vanz].Color->type       = cf::D;     
	    v[vanz].Color->SetParticleArg(1,2);     
	    v[vanz].Color->SetStringArg('1','2');     
	    
	    v[vanz].nlf     = 1;
	    v[vanz].Lorentz = new Lorentz_Function; 
	    
	    v[vanz].Lorentz->type = lf::VVSS;     
	    v[vanz].Lorentz->SetParticleArg(0,3);     
	    
	    v[vanz].on      = 1;
	    vanz++;
	  }
	  
	  if (flavZ.IsOn()) {
	    // W - D - U - Z  
	    
	    v[vanz].in[0] = flavW;
	    v[vanz].in[1] = flav1.Bar();
	    v[vanz].in[2] = flav2;
	    v[vanz].in[3] = flavZ;
	    
	    v[vanz].nleg     = 4;
	    
	    kcpl0 = num_2;
	    kcpl1 = kcpl0;
	    
	    v[vanz].cpl[0]  = kcpl0.Value();
	    v[vanz].cpl[1]  = kcpl1.Value();
	    v[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	    v[vanz].cpl[2]  = 0.;v[vanz].cpl[3]  = 0.;
	    
	    v[vanz].ncf   = 1;
	    v[vanz].Color = new Color_Function; 
	    
	    v[vanz].Color->type       = cf::D;     
	    v[vanz].Color->SetParticleArg(1,2);     
	    v[vanz].Color->SetStringArg('1','2');     
	    
	    v[vanz].nlf     = 1;
	    v[vanz].Lorentz = new Lorentz_Function; 
	    
	    v[vanz].Lorentz->type = lf::VVSS;     
	    v[vanz].Lorentz->SetParticleArg(0,3);     
	    
	    v[vanz].on      = 1;
	    vanz++;
	  }
	}
      }
    }
  }
}

Kabbala Model_SQCD::K_Z_D(short int i,short int j)       
{   
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("Z^{")+string(hi)+string(hj)+string("}_D"),SpsD.Zd(i,j));
}  
Kabbala Model_SQCD::K_Z_U(short int i,short int j)       
{   
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("Z^{")+string(hi)+string(hj)+string("}_U"),SpsU.Zu(i,j));
}  
Kabbala Model_SQCD::K_w_S(short int i,short int j)
{   
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("w^{")+string(hi)+string(hj)+string("}_S"),SpsU.w_S(i,j));
}  
Kabbala Model_SQCD::K_u_S(short int i,short int j)
{   
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("u^{")+string(hi)+string(hj)+string("}_S"),SpsU.u_S(i,j));
}  
Kabbala Model_SQCD::K_e_S(short int i,short int j)
{   
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("e^{")+string(hi)+string(hj)+string("}_S"),SpsD.e_S(i,j));
}  
Kabbala Model_SQCD::K_d_S(short int i,short int j)
{   
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("d^{")+string(hi)+string(hj)+string("}_S"),SpsD.d_S(i,j));
}  
int Model_SQCD::gen_sUp(Flavour fl)
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
int Model_SQCD::gen_sDown(Flavour fl)
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


