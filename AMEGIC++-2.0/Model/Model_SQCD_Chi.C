#include "Model_SQCD_Chi.H"
#include <stdio.h>

using namespace AMEGIC;
using namespace APHYTOOLS;
using namespace std;

Model_SQCD_Chi::~Model_SQCD_Chi()
{
  if (ini_sm==1) delete mosm;
  if (ini_ch==1) delete moch;
  if (ini_sq==1) delete mosq;
}

void Model_SQCD_Chi::Init()
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
  if (mosq==NULL) {
    //SQCD part
    mosq = new Model_SQCD;
    mosq->Init();
    ini_sq = 1;
  }

  g1     = Kabbala(string("e"),sqrt(4.*M_PI*mosm->Aqed()));
  g2     = Kabbala(string("e/sin\\theta_W"), g1.Value()/mosm->SinTW());
  PL     = Kabbala(string("P_L"),1.);
  PR     = Kabbala(string("P_R"),1.);
  M_I    = Kabbala(string("i"),Complex(0.,1.));
  root2  = Kabbala(string("\\sqrt{2}"),sqrt(2.));
  K_zero = Kabbala(string("zero"),0.);
  num_2 = Kabbala(string("2"),2.);    	
  num_3 = Kabbala(string("3"),3.);    	
}

void Model_SQCD_Chi::Init(Model_SM* _mosm,Model_Chi* _moch,Model_SQCD* _mosq)
{
  mosm = _mosm;
  moch = _moch;
  mosq = _mosq;

  Init();
}

void Model_SQCD_Chi::c_FFS(Single_Vertex* v,int& vanz){
  
  Kabbala kcpl0,kcpl1; 
  
  //quark - squark - neutralino
  
  for (short int j=43;j<47;j++) {
    Flavour flneu = Flavour(kf::code(j));
    if (flneu.ison()) {
      //uptypes 
      for (short int k=2;k<7;k+=2) {
	Flavour flav1 = Flavour(kf::code(k));
	for (short int i=51;i<57;i++) {
	  Flavour flav2 = Flavour(kf::code(i));
	  if (flav1.ison() && flav2.ison()) {
	    v[vanz].in[0] = flav1;
	    v[vanz].in[1] = flav2;
	    v[vanz].in[2] = flneu;
	    
	    /*Kabbala K_uI = Kabbala(string("u^I"),
	      flav1.yuk()*sqrt(2.)/(mosq->K_v2()).Value());*/

	    Kabbala K_u1 = Kabbala(string("u^1"),
				   Flavour(kf::code(2)).yuk()*sqrt(2.)/(mosq->K_v2()).Value());
	    Kabbala K_u2 = Kabbala(string("u^2"),
				   Flavour(kf::code(4)).yuk()*sqrt(2.)/(mosq->K_v2()).Value());
	    Kabbala K_u3 = Kabbala(string("u^3"),
				   Flavour(kf::code(6)).yuk()*sqrt(2.)/(mosq->K_v2()).Value());
	    
	    kcpl0 = M_I*(((g1*root2*num_2)/
			  (K_cosTW()*num_3))*K_Z_U((k-2)/2+3,i-51)*K_Z_N(0,j-43)-
			 (K_u1*K_Z_U(0,i-51)+
			  K_u2*K_Z_U(1,i-51)+
			  K_u3*K_Z_U(2,i-51))*K_Z_N(3,j-43));
	    
	    kcpl1 = M_I*(-g2/(K_cosTW()*root2)*K_Z_U((k-2)/2,i-51)*
			 (K_Z_N(0,j-43)*(K_sinTW()/num_3)+K_Z_N(1,j-43)*K_cosTW())-
		       (K_u1*K_Z_U(3,i-51)+
			K_u2*K_Z_U(4,i-51)+
			K_u3*K_Z_U(5,i-51))*K_Z_N(3,j-43));
	    
	    v[vanz].cpl[0] = kcpl0.Value();
	    v[vanz].cpl[1] = kcpl1.Value();
	    v[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();
	    v[vanz].cpl[2]  = 0.;v[vanz].cpl[3]  = 0.;

	    v[vanz].ncf   = 1;
	    v[vanz].Color = new Color_Function; 
	    
	    v[vanz].Color->type       = cf::D;     
	    v[vanz].Color->SetParticleArg(0,1);     
	    v[vanz].Color->SetStringArg('0','1');     
	  
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

  for (short int j=43;j<47;j++) {
    Flavour flneu = Flavour(kf::code(j));
    if (flneu.ison()) {
      //downtypes 
      for (short int k=1;k<6;k+=2) {
	  Flavour flav1 = Flavour(kf::code(k));
	  for (short int i=61;i<67;i++) {
	    Flavour flav2 = Flavour(kf::code(i));
	    if (flav1.ison() && flav2.ison()) {
	      v[vanz].in[0] = flav1;
	      v[vanz].in[1] = flav2;
	      v[vanz].in[2] = flneu;

	      Kabbala K_dI = Kabbala(string("d^I"),
		-flav1.yuk()*sqrt(2.)/(mosq->K_v1()).Value());
	      Kabbala K_d1 = Kabbala(string("d^1"),
				     -Flavour(kf::code(1)).yuk()*sqrt(2.)/(mosq->K_v1()).Value());
	      Kabbala K_d2 = Kabbala(string("d^2"),
				     -Flavour(kf::code(3)).yuk()*sqrt(2.)/(mosq->K_v1()).Value());
	      Kabbala K_d3 = Kabbala(string("d^3"),
				     -Flavour(kf::code(5)).yuk()*sqrt(2.)/(mosq->K_v1()).Value());
				     
	      kcpl0 = M_I*((-(g1*root2)/
			    (K_cosTW()*num_3))*K_Z_D((k-1)/2+3,i-61)*K_Z_N(0,j-43)+
			   (K_d1*K_Z_D(0,i-61)+
			   K_d2*K_Z_D(1,i-61)+
			    K_d3*K_Z_D(2,i-61))/*(K_dI*K_Z_D((k-1)/2,i-61))*/*K_Z_N(2,j-43));
	      
	      kcpl1 = M_I*(-g2/(K_cosTW()*root2)*K_Z_D((k-1)/2,i-61)*
			   (K_Z_N(0,j-43)*(K_sinTW()/num_3)-K_Z_N(1,j-43)*K_cosTW())+
			   (K_d1*K_Z_D(3,i-61)+
			   K_d2*K_Z_D(4,i-61)+
			    K_d3*K_Z_D(5,i-61))/*(K_dI*K_Z_D((k-1)/2+3,i-61))*/*K_Z_N(2,j-43));
	      
	      v[vanz].cpl[0] = kcpl0.Value();
	      v[vanz].cpl[1] = kcpl1.Value();
	      v[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();
	      v[vanz].cpl[2] = 0.;v[vanz].cpl[3]  = 0.;

	      v[vanz].ncf   = 1;
	      v[vanz].Color = new Color_Function; 
	      
	      v[vanz].Color->type       = cf::D;     
	      v[vanz].Color->SetParticleArg(0,1);     
	      v[vanz].Color->SetStringArg('0','1');     
	      
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

  //d-quark - Chargino - sup
  for (short int i=1;i<6;i+=2) {
    Flavour flav1 = Flavour(kf::code(i));
      
    Kabbala K_dI = Kabbala(string("d^I"),
			   -flav1.yuk()/(mosq->K_v1()).Value()*sqrt(2.));
    
    for (short int j=41;j<43;j++) {
      Flavour flav2 = Flavour(kf::code(j));
      for (short int k=51;k<57;k++) {
	Flavour flav3 = Flavour(kf::code(k));
	if (flav1.ison() && flav2.ison() && flav3.ison()) {
	  v[vanz].in[0] = flav1;
	  v[vanz].in[1] = flav3;
	  v[vanz].in[2] = flav2;
	
	  Kabbala K_uI = Kabbala(string("u^I"),Flavour(kf::code(2*gen_sUp(flav3)+2)).yuk()*
				 sqrt(2.)/(mosq->K_v2()).Value());

	  kcpl0 = -M_I*K_dI*K_Z_MI(1,j-41)*K_Z_U(gen_sUp(flav3),k-51)*
	    K_CKM((i-1)/2,gen_sUp(flav3));
	  
	  kcpl1 = M_I*(-g2*K_Z_PL(0,j-41)*K_Z_U(gen_sUp(flav3),k-51)+
		       K_uI*K_Z_PL(1,j-41)*K_Z_U(gen_sUp(flav3)+3,k-51))*
		       K_CKM((i-1)/2,gen_sUp(flav3));

	  v[vanz].cpl[0] = kcpl0.Value();
	  v[vanz].cpl[1] = kcpl1.Value();
	  v[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();
	  v[vanz].cpl[2]  = 0.;v[vanz].cpl[3]  = 0.;
	  
	  v[vanz].ncf   = 1;
	  v[vanz].Color = new Color_Function; 
	
	  v[vanz].Color->type       = cf::D;     
	  v[vanz].Color->SetParticleArg(0,1);     
	  v[vanz].Color->SetStringArg('0','1');     
	  
	  v[vanz].nlf     = 1;
	  v[vanz].Lorentz = new Lorentz_Function; 
	      
	  v[vanz].Lorentz->type       = lf::FFS;     
	  
	  v[vanz].on      = 1;
	  vanz++;
	}
     }
    }
  }

  //u-quark - Chargino - sdown
  for (short int i=2;i<7;i+=2) {
    Flavour flav1 = Flavour(kf::code(i));
    
    Kabbala K_uJ = Kabbala(string("u^J"),flav1.yuk()*sqrt(2.)/(mosq->K_v2()).Value());
    
    for (short int j=41;j<43;j++) {
      Flavour flav2 = Flavour(kf::code(j));
      for (short int k=61;k<67;k++) {
	Flavour flav3 = Flavour(kf::code(k));
	if (flav1.ison() && flav2.ison() && flav3.ison()) {
	  v[vanz].in[0] = flav1;
	  v[vanz].in[1] = flav3;
	  v[vanz].in[2] = flav2.bar();

	  Kabbala K_dI = Kabbala(string("d^I"),
				 -Flavour(kf::code(2*gen_sDown(flav3)+1)).yuk()*sqrt(2.)/
				 (mosq->K_v1()).Value());
	  
	  kcpl0 = M_I*K_uJ*K_Z_D(gen_sDown(flav3),k-61)*K_Z_PL(1,j-41)*
	    K_CKM_conj(gen_sDown(flav3),(i-2)/2);
	  
	  kcpl1 = -M_I*(g2*K_Z_D(gen_sDown(flav3),k-61)*K_Z_MI(0,j-41)+
			K_dI*K_Z_D(gen_sDown(flav3)+3,k-61)*K_Z_MI(1,j-41))*
			K_CKM_conj(gen_sDown(flav3),(i-2)/2);
	  
	  v[vanz].cpl[0] = kcpl0.Value();
	  v[vanz].cpl[1] = kcpl1.Value();
	  v[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();
	  v[vanz].cpl[2] = 0.;v[vanz].cpl[3]  = 0.;
	  
	  v[vanz].ncf   = 1;
	  v[vanz].Color = new Color_Function; 
	  
	  v[vanz].Color->type       = cf::D;     
	  v[vanz].Color->SetParticleArg(0,1);     
	  v[vanz].Color->SetStringArg('0','1');     
	  
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

inline Kabbala Model_SQCD_Chi::K_CKM_conj(short int i,short int j)       
{   
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("V_{")+string(hi)+string(hj)+string("}^\\ti"),
		 conj((mosm->K_CKM(i,j)).Value()));
} 



















