#include <stdio.h>
#include "Model_Higgs.H"
//#include "Run_Parameter.H"
#include "MathTools.H"

using namespace AMEGIC;
using namespace APHYTOOLS;
using namespace std;

// Missing : Running Yukawa couplings !!!



Model_Higgs::~Model_Higgs()
{
  if (ini_sm==1)    delete mosm;
}

void Model_Higgs::Init()
{
  
  if (mosm==NULL) {
   //SM part
    mosm = new Model_SM;
    mosm->Init();
    ini_sm = 1;
  }
  if (isa!=NULL) SpHiggs.Interface(isa);
 else
    SpHiggs.Init();
 
 
  g1    = Kabbala(string("e"),sqrt(4.*M_PI*mosm->Aqed()));
  g2    = Kabbala(string("e/sin\\theta_W"), g1.Value()/mosm->SinTW());
  PL     = Kabbala(string("P_L"),1.);
  PR     = Kabbala(string("P_R"),1.);
  M_I    = Kabbala(string("i"),Complex(0.,1.));
  root2  = Kabbala(string("\\sqrt{2}"),sqrt(2.));
  K_zero = Kabbala(string("zero"),0.);
  num_2 = Kabbala(string("2"),2.);    	
  num_3 = Kabbala(string("3"),3.);    	
  num_4 = Kabbala(string("4"),4.);    		
}

void Model_Higgs::Init(Model_SM* _mosm,Isajet* _isa)
{
  mosm    = _mosm;
  isa = _isa;
  Init();
}

void Model_Higgs::c_FFS(Single_Vertex* v,int& vanz)
{
  Flavour flHmin(kf::Hmin);
  Flavour flH0(kf::H0);
  Flavour flh0(kf::h0);
  Flavour flA0(kf::A0);
  
  // l(e)/q(d) -> h0/H0 + l(e)/q(d)
  
  for(short int i=31;i<33;i++) {
    
    Flavour flav = Flavour(kf::code(i));
    Kabbala kcpl0,kcpl1;

    if(flav.ison()) {
      
      for(short int j=1;j<17;j++) {
	if (j==7) j=11;
	Flavour fl1 = Flavour(kf::code(j)); 
       	if(fl1.ison() && fl1.isfermion() && fl1.isdowntype() && (fl1.yuk() > 1.)) {
	  
	  v[vanz].in[0] = fl1; 
	  v[vanz].in[1] = flav;
	  v[vanz].in[2] = fl1; 
	  
	  kcpl0 = -M_I/K_v1()*K_yuk(fl1)*K_Z_R(0,i-31);
	  kcpl1 = kcpl0;

	  v[vanz].cpl[0]  = kcpl0.Value();
	  v[vanz].cpl[1]  = kcpl1.Value();
	  v[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	  v[vanz].cpl[2]  = 0.;v[vanz].cpl[3]  = 0.;
	  v[vanz].on      = 1;

	  v[vanz].ncf   = 1;
	  v[vanz].Color = new Color_Function; 
	  
	  if (fl1.strong()) {
	    v[vanz].Color->type       = cf::D;     
	    v[vanz].Color->SetParticleArg(0,2);     
	    v[vanz].Color->SetStringArg('0','2');     
	  }
	  else v[vanz].Color->type = cf::None; 
	
	  v[vanz].nlf     = 1;
	  v[vanz].Lorentz = new Lorentz_Function; 
	  
	  v[vanz].Lorentz->type       = lf::FFS;     
	  	  
	  vanz++;
	  //checked FK & RK	  
	}
      }
    }
  
    // q(u) -> h0/H0 + q(u)
    
    for(short int t=1;t<17;t++) {
      if (t==7) t=11;
      Flavour fl1 = Flavour(kf::code(t)); 
      if(fl1.ison() && fl1.isquark() && fl1.isuptype() && (fl1.yuk() > 1.)) {
	
	v[vanz].in[0] = fl1; 
	v[vanz].in[1] = flav;
	v[vanz].in[2] = fl1; 

	kcpl0 = -M_I/K_v2()*K_yuk(fl1)*K_Z_R(1,i-31);
	kcpl1 = kcpl0;

	v[vanz].cpl[0]  = kcpl0.Value();
	v[vanz].cpl[1]  = kcpl1.Value();
	v[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	v[vanz].cpl[2]  = 0.;v[vanz].cpl[3]  = 0.;
	v[vanz].on      = 1;
	
	v[vanz].ncf   = 1;
	v[vanz].Color = new Color_Function; 
	  
	v[vanz].Color->type       = cf::D;     
	v[vanz].Color->SetParticleArg(0,2);     
	v[vanz].Color->SetStringArg('0','2');     
	
	v[vanz].nlf     = 1;
	v[vanz].Lorentz = new Lorentz_Function; 
	
	v[vanz].Lorentz->type       = lf::FFS;     
	v[vanz].Lorentz->SetParticleArg(1);     
	
	vanz++;
	//checked FK & RK
      }
    }
  }
  
  // l(e)/q(d) -> A0 + l(e)/q(d)
  
  if(flA0.ison()) {

    Kabbala kcpl0,kcpl1;
 
    for(short int k=1;k<17;k++) {
      if (k==7) k=11;
      Flavour fl1 = Flavour(kf::code(k)); 
      if(fl1.ison() && fl1.isfermion() && fl1.isdowntype() && (fl1.yuk() > 1.)) {
		
	v[vanz].in[0] = fl1; 
	v[vanz].in[1] = flA0;
	v[vanz].in[2] = fl1; 


	Kabbala K_rpa_Mass;
	K_rpa_Mass = Kabbala(string("M_{")+fl1.texname()+string("}"),fl1.yuk());
	//	K_rpa_Mass = Kabbala(string("M_{")+fl1.texname()+string("}"),
	//			     rpa.consts.Mass(fl1,sqr(flA0.mass())));

	kcpl0 = -K_rpa_Mass/K_v1()*K_Z_H(0,0);
	kcpl1 = -kcpl0;

	v[vanz].cpl[0]  = kcpl0.Value();
	v[vanz].cpl[1]  = kcpl1.Value();
	v[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	v[vanz].cpl[2]  = 0.;v[vanz].cpl[3]  = 0.;
	v[vanz].on      = 1;
	
	v[vanz].ncf   = 1;
	v[vanz].Color = new Color_Function; 
	  
	if (fl1.strong()) {
	  v[vanz].Color->type       = cf::D;     
	  v[vanz].Color->SetParticleArg(0,2);     
	  v[vanz].Color->SetStringArg('0','2');     
	  }
	else v[vanz].Color->type = cf::None; 
	
	v[vanz].nlf     = 1;
	v[vanz].Lorentz = new Lorentz_Function; 
	  
	v[vanz].Lorentz->type       = lf::FFS;     
	v[vanz].Lorentz->SetParticleArg(1);     
	
	vanz++;
	//checked FK & RK	
      }
    }
    
    // q(u) -> A0 + q(u)
    
    for(short int z=1;z<7;z++) { 
      Flavour fl1 = Flavour(kf::code(z)); 
      if(fl1.ison() && fl1.isquark() && fl1.isuptype() && (fl1.yuk() > 1.)) {
	
	v[vanz].in[0] = fl1; 
	v[vanz].in[1] = flA0;
	v[vanz].in[2] = fl1; 

	kcpl0 = -K_yuk(fl1)/K_v2()*K_Z_H(1,0);
	kcpl1 = -kcpl0;

	v[vanz].cpl[0]  = kcpl0.Value();
	v[vanz].cpl[1]  = kcpl1.Value();
	v[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	v[vanz].cpl[2]  = 0.;v[vanz].cpl[3]  = 0.;
	v[vanz].on      = 1;
	
	v[vanz].ncf   = 1;
	v[vanz].Color = new Color_Function; 
	
	v[vanz].Color->type       = cf::D;     
	v[vanz].Color->SetParticleArg(0,2);     
	v[vanz].Color->SetStringArg('0','2');     
      
	v[vanz].nlf     = 1;
	v[vanz].Lorentz = new Lorentz_Function; 
	
	v[vanz].Lorentz->type       = lf::FFS;     
	v[vanz].Lorentz->SetParticleArg(1);     

	vanz++;
	//checked FK & RK
      }
    }
  } 
  
  if(flHmin.ison()) {
    
    // l(e) -> H- + nu(e)
    
    Kabbala kcpl0,kcpl1;

    for(short int i=11;i<17;i++) {
      Flavour fl1 = Flavour(kf::code(i)); 
      if(fl1.ison() && fl1.islepton() && fl1.isdowntype() && (fl1.yuk() > 1.)) {
	Flavour fl2 = Flavour(kf::code(i=i+1));
	if(fl2.ison() && fl2.islepton() && fl2.isuptype() ) {
	  
	  v[vanz].in[0] = fl1;   
	  v[vanz].in[1] = flHmin;
	  v[vanz].in[2] = fl2;   

	  kcpl1 = M_I/K_v1()*root2*K_yuk(fl1)*K_Z_H(0,0);
	  kcpl0 = K_zero;	 	  

	  v[vanz].cpl[0]  = kcpl0.Value();
 	  v[vanz].cpl[1]  = kcpl1.Value();
	  v[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();	  
	  v[vanz].cpl[2]  = 0.;v[vanz].cpl[3]  = 0.;
	  v[vanz].on      = 1;

	  v[vanz].ncf   = 1;
	  v[vanz].Color = new Color_Function; 
	  	  
	  v[vanz].Color->type       = cf::None;     
	  
	  v[vanz].nlf     = 1;
	  v[vanz].Lorentz = new Lorentz_Function; 
	  
	  v[vanz].Lorentz->type       = lf::FFS;     
	  v[vanz].Lorentz->SetParticleArg(1);     
	  
	  vanz++;
	  //checked FK & RK	  
	}
      }
    }  
    
    // q(d) -> H- + q(u) plus Generation Mixing !
    
    for(short int j=1;j<7;j++) {
      Flavour fl1 = Flavour(kf::code(j));
      Kabbala kcpl0,kcpl1;
      if(fl1.ison() && fl1.isquark() && fl1.isdowntype()) {
	for(short int k=1;k<7;k++) {
	  Flavour fl2 = Flavour(kf::code(k));
	  if(fl2.ison() && fl2.isquark() && fl2.isuptype() && 
	     ((fl1.yuk() > 1.) || (fl2.yuk() > 1.)) ) {
	    
            int geni=(fl1.kfcode()-1)/2; //downtype
	    int genj=(fl2.kfcode()-2)/2; //uptype
	   
	    v[vanz].in[0] = fl1;
	    v[vanz].in[1] = flHmin;
	    v[vanz].in[2] = fl2;
	   
	    //unbedingt noch einmal checken !!!!!!!!!!!!!!!

	    kcpl0 = M_I/K_v2()*root2*K_yuk(fl2)*K_Z_H(1,0)*K_CKM(geni,genj);
	    kcpl1 = M_I/K_v1()*root2*K_yuk(fl1)*K_Z_H(0,0)*K_CKM(geni,genj);
	    	   
	    v[vanz].cpl[0]  = kcpl0.Value();
	    v[vanz].cpl[1]  = kcpl1.Value();
	    v[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();	  	    
	    v[vanz].cpl[2]  = 0.;v[vanz].cpl[3]  = 0.;
	    v[vanz].on      = 1;

	    v[vanz].ncf   = 1;
	    v[vanz].Color = new Color_Function; 
	  	    
	    v[vanz].Color->type       = cf::D;     
	    v[vanz].Color->SetParticleArg(0,2);     
	    v[vanz].Color->SetStringArg('0','2');     
	    
	    v[vanz].nlf     = 1;
	    v[vanz].Lorentz = new Lorentz_Function; 
	    
	    v[vanz].Lorentz->type       = lf::FFS;     
	    v[vanz].Lorentz->SetParticleArg(1);     
	    
	    vanz++;
	    //checked FK & RK
	  }
	}
      }  
    }      
  }
}

void Model_Higgs::c_VVS(Single_Vertex* v,int& vanz) 
{
    
  Flavour flW(kf::W);
  Flavour flZ(kf::Z);
  Flavour flHmin(kf::Hmin);    
  Flavour flh0(kf::h0);
  Flavour flH0(kf::H0);
  Flavour flA0(kf::A0);
  Flavour flPhoton(kf::photon);
  
  Kabbala kcpl0,kcpl1;
  
  //W- -> h0,H0 + W-

  if(flW.ison()) {
    
    for(short int i=31; i<33;i++){
      Flavour flav = Flavour(kf::code(i));
      if(flav.ison()) {
	
	v[vanz].in[0] = flW;
	v[vanz].in[1] = flav;
	v[vanz].in[2] = flW;

	kcpl0 = M_I/num_2*g2*g2*(K_v1()*K_Z_R(0,i-31)+K_v2()*K_Z_R(1,i-31));
	kcpl1 = kcpl0;

	v[vanz].cpl[0]  = kcpl0.Value();
	v[vanz].cpl[1]  = kcpl1.Value();
	v[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();	  	    	
	v[vanz].cpl[2]  = 0.;v[vanz].cpl[3]  = 0.;
	v[vanz].on      = 1;
	
	v[vanz].ncf   = 1;
	v[vanz].Color = new Color_Function; 
	
	v[vanz].Color->type       = cf::None;     
	
	v[vanz].nlf     = 1;
	v[vanz].Lorentz = new Lorentz_Function; 
	
	v[vanz].Lorentz->type = lf::Gab;     
	v[vanz].Lorentz->SetParticleArg(0,2);     
	
	vanz++;
	//checked FK & RK
      }
    }
  }
  
  //  Z -> h0,H0 + Z  
  
  if(flZ.ison()){
    for(short int i=31; i<33;i++){
      Flavour flav = Flavour(kf::code(i));
      if(flav.ison()) {
	
 	v[vanz].in[0] = flZ;
	v[vanz].in[1] = flav;
	v[vanz].in[2] = flZ;
	
	kcpl0 = M_I/(K_cosTW()*K_cosTW()*num_2)*g2*g2*
	  (K_v1()*K_Z_R(0,i-31)+K_v2()*K_Z_R(1,i-31));

	kcpl1 = kcpl0;

	v[vanz].cpl[0]  = kcpl0.Value();
	v[vanz].cpl[1]  = kcpl1.Value();
	v[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();	  	    	
	v[vanz].cpl[2]  = 0.;v[vanz].cpl[3]  = 0.;
	v[vanz].on      = 1;
	
	v[vanz].ncf   = 1;
	v[vanz].Color = new Color_Function; 
	
	v[vanz].Color->type       = cf::None;     
	
	v[vanz].nlf     = 1;
	v[vanz].Lorentz = new Lorentz_Function; 
	
	v[vanz].Lorentz->type = lf::Gab;     
	v[vanz].Lorentz->SetParticleArg(0,2);     
	
	vanz++;
	//checked FK & RK & SS
      }
    }  
  }
}


void Model_Higgs::c_SSS(Single_Vertex* v,int& vanz) 
{
  Flavour flHmin(kf::Hmin);    
  Flavour flh0(kf::h0);
  Flavour flH0(kf::H0);
  Flavour flA0(kf::A0);
  Kabbala kcpl0,kcpl1;

  // h0,H0 -> H+ + H-
  
  if(flHmin.ison()) {  
    for(short int i=31;i<33;i++) {
      Flavour flav = Flavour(kf::code(i)); 
      if(flav.ison()) {
	
	v[vanz].in[0] = flHmin;
	v[vanz].in[1] = flav;
	v[vanz].in[2] = flHmin;
	
	//unbedingt noch einmal checken !!!!

	kcpl0 = -M_I*g2*g2*(K_A_H(0,0)*K_B_R(i-31)/(K_cosTW()*K_cosTW()*num_4)+
			    K_v0()*K_A_P(i-31,0)/num_2);
	kcpl1 = kcpl0;
	
	v[vanz].cpl[0]  = kcpl0.Value();
	v[vanz].cpl[1]  = kcpl1.Value();
	v[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();	  	    		
	v[vanz].cpl[2]  = 0.;v[vanz].cpl[3]  = 0.;
	
	v[vanz].ncf   = 1;
	v[vanz].Color = new Color_Function; 
	
	v[vanz].Color->type = cf::None; 
	
	v[vanz].nlf     = 1;
	v[vanz].Lorentz = new Lorentz_Function; 
	
	v[vanz].Lorentz->type       = lf::SSS;     
	
	v[vanz].on      = 1;
	vanz++;
	//checked FK & RK
      }  
      
    }
  }
  
  // h0 -> h0/H0 + h0/H0 

    if(flh0.ison()) { 
      for(short int j=31;j<33;j++) {
	Flavour flav2 = Flavour(kf::code(j));
	if(flav2.ison()) {
	  for(short int k=31;k<33;k++) {
	    Flavour flav3 = Flavour(kf::code(k));
	    if(flav3.ison()) {
	      if(flav2 != flH0 || flav3 != flh0){

	      v[vanz].in[0] = flh0;
	      v[vanz].in[1] = flav2;      
	      v[vanz].in[2] = flav3;
	      
	      Kabbala C1 = K_Z_R(0,0)*K_Z_R(0,j-31)*K_Z_R(0,k-31);
	      Kabbala C2 = K_Z_R(1,0)*K_Z_R(1,j-31)*K_Z_R(1,k-31);
	      Kabbala C3 = K_Z_R(0,0)*K_Z_R(0,j-31)*K_Z_R(1,k-31);
	      Kabbala C4 = K_Z_R(0,0)*K_Z_R(1,j-31)*K_Z_R(0,k-31);
	      Kabbala C5 = K_Z_R(1,0)*K_Z_R(0,j-31)*K_Z_R(0,k-31);
	      Kabbala C6 = K_Z_R(1,0)*K_Z_R(1,j-31)*K_Z_R(0,k-31);
	      Kabbala C7 = K_Z_R(1,0)*K_Z_R(0,j-31)*K_Z_R(1,k-31);
	      Kabbala C8 = K_Z_R(0,0)*K_Z_R(1,j-31)*K_Z_R(1,k-31);
	      
	      kcpl0 = -M_I*g2*g2/(K_cosTW()*K_cosTW()*num_4)*
		(K_v1()*(C1*num_3-C6-C7-C8)+K_v2()*(C2*num_3-C3-C4-C5));
	      kcpl1 = kcpl0;

	      v[vanz].cpl[0]  = kcpl0.Value();
	      v[vanz].cpl[1]  = kcpl1.Value();
	      v[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();	  	    		
	      v[vanz].cpl[2]  = 0.;v[vanz].cpl[3]  = 0.;
	
	      v[vanz].ncf   = 1;
	      v[vanz].Color = new Color_Function; 
	      
	      v[vanz].Color->type = cf::None; 
    
	      v[vanz].nlf     = 1;
	      v[vanz].Lorentz = new Lorentz_Function; 
	      
	      v[vanz].Lorentz->type       = lf::SSS;     
	      
	      v[vanz].on      = 1;
	      vanz++;
	      //checked FK & RK
	      }
	    }
	  }
	} 
      }
    }
    
  
  //H0 -> H0 + H0

    if(flH0.ison()) { 
      v[vanz].in[0] = flH0;
      v[vanz].in[1] = flH0;      
      v[vanz].in[2] = flH0;
      
      Kabbala C1 = K_Z_R(0,1)*K_Z_R(0,1)*K_Z_R(0,1);
      Kabbala C2 = K_Z_R(1,1)*K_Z_R(1,1)*K_Z_R(1,1);
      Kabbala C3 = K_Z_R(1,1)*K_Z_R(0,1)*K_Z_R(0,1);
      Kabbala C4 = K_Z_R(0,1)*K_Z_R(1,1)*K_Z_R(1,1);
      
      //noch einmal checken !!

      kcpl0 = -M_I*g2*g2/(K_cosTW()*K_cosTW()*num_4)*
	((C1-C4)*K_v1()*num_3+(C2-C3)*K_v2()*num_3);
      kcpl1 = kcpl0;

      v[vanz].cpl[0]  = kcpl0.Value();
      v[vanz].cpl[1]  = kcpl1.Value();
      v[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();	  	    		
      v[vanz].cpl[2]  = 0.;v[vanz].cpl[3]  = 0.;
      
      v[vanz].ncf   = 1;
      v[vanz].Color = new Color_Function; 
      
      v[vanz].Color->type = cf::None; 
      
      v[vanz].nlf     = 1;
      v[vanz].Lorentz = new Lorentz_Function; 
      
      v[vanz].Lorentz->type       = lf::SSS;     
      
      v[vanz].on      = 1;
      vanz++;
    }
      
  // A0 -> h0/H0 + A0
    
    if(flA0.ison()) {
      for(short int i=31;i<33;i++) {
	Flavour flav = Flavour(kf::code(i));
	if(flav.ison()) {
	  
	v[vanz].in[0] = flA0;
	v[vanz].in[1] = flav;
	v[vanz].in[2] = flA0;
	
	//noch einmal checken !!!

	kcpl0 = -M_I/(K_cosTW()*K_cosTW()*num_4)*g2*g2*K_A_H(0,0)*K_B_R(i-31);
	kcpl1 = kcpl0;
	
	v[vanz].cpl[0]  = kcpl0.Value();
	v[vanz].cpl[1]  = kcpl1.Value();
	v[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();	  	    		
	v[vanz].cpl[2]  = 0.;v[vanz].cpl[3]  = 0.;
	
	v[vanz].ncf   = 1;
	v[vanz].Color = new Color_Function; 
	
	v[vanz].Color->type = cf::None; 
	
	v[vanz].nlf     = 1;
	v[vanz].Lorentz = new Lorentz_Function; 
	
	v[vanz].Lorentz->type       = lf::SSS;     
	
	v[vanz].on      = 1;
	vanz++;
	//checked FK & RK	
      }
    }
  }
 
}


void Model_Higgs::c_SSV(Single_Vertex* v,int& vanz)
{   
  Flavour flW(kf::W);
  Flavour flZ(kf::Z);
  Flavour flHmin(kf::Hmin);
  Flavour flh0(kf::h0);
  Flavour flH0(kf::H0);
  Flavour flA0(kf::A0);
  Flavour flPhoton(kf::photon);
  Kabbala kcpl0,kcpl1;

  // A0 -> Z + h0/H0 
  
  if(flZ.ison()) {
    for(short int i=31; i<33;i++) {
      Flavour flav = Flavour(kf::code(i)); 
      if(flav.ison() && flA0.ison()) {
	v[vanz].in[0] = flA0;
	v[vanz].in[1] = flZ;
	v[vanz].in[2] = flav;

	kcpl0 = g2/(K_cosTW()*num_2)*K_A_M(i-31,0);
	kcpl1 = kcpl0;



	v[vanz].cpl[0]  = kcpl0.Value();
	v[vanz].cpl[1]  = kcpl1.Value();
	v[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();	  	    		
	v[vanz].cpl[2]  = 0.;v[vanz].cpl[3]  = 0.;
	
	v[vanz].ncf   = 1;
	v[vanz].Color = new Color_Function; 
	
	v[vanz].Color->type = cf::None; 
 
	v[vanz].nlf     = 1;
	v[vanz].Lorentz = new Lorentz_Function; 
      
	v[vanz].Lorentz->type       = lf::SSV;     
	v[vanz].Lorentz->SetParticleArg(0,2,1);     	

	v[vanz ].on      = 1;
	vanz++;
	//checked FK & RK
      }  
    }
    
    // H- -> Z + H-
    
    if(flHmin.ison()) {
      v[vanz].in[0] = flHmin;
      v[vanz].in[1] = flZ;
      v[vanz].in[2] = flHmin;
     

      Kabbala cot2TW = Kabbala(string("cot2\\theta_W"),
			       (CosTW()*CosTW()-SinTW()*SinTW())/(2.*SinTW()*CosTW()));

      kcpl0 = M_I*g1*cot2TW;
      kcpl1 = kcpl0;

      v[vanz].cpl[0]  = kcpl0.Value();
      v[vanz].cpl[1]  = kcpl1.Value();
      v[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();	  	    		
      v[vanz].cpl[2]  = 0.;v[vanz].cpl[3]  = 0.;
      
      v[vanz].ncf   = 1;
      v[vanz].Color = new Color_Function; 
	        
      v[vanz].Color->type = cf::None; 
      
      v[vanz].nlf     = 1;
      v[vanz].Lorentz = new Lorentz_Function; 
      
      v[vanz].Lorentz->type       = lf::SSV;     
      v[vanz].Lorentz->SetParticleArg(0,2,1);     	
      
      v[vanz].on      = 1;
      vanz++;
      //checked FK & RK
    }    
  }
  
  // H- -> Photon + H- 
  
  if(flPhoton.ison() && flHmin.ison()) {
    
    v[vanz].in[0] = flHmin;
    v[vanz].in[1] = flPhoton;
    v[vanz].in[2] = flHmin;
   
    //sollte so stimmen !

    kcpl0 = -M_I*g1;
    kcpl1 = kcpl0;

    v[vanz].cpl[0]  = kcpl0.Value();
    v[vanz].cpl[1]  = kcpl1.Value();
    v[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();	  	    		
    v[vanz].cpl[2]  = 0.;v[vanz].cpl[3]  = 0.;
    
    v[vanz].ncf   = 1;
    v[vanz].Color = new Color_Function; 
    
    v[vanz].Color->type = cf::None; 
    
    v[vanz].nlf     = 1;
    v[vanz].Lorentz = new Lorentz_Function; 
    
    v[vanz].Lorentz->type       = lf::SSV;     
    v[vanz].Lorentz->SetParticleArg(0,2,1);     	
          
    v[vanz].on      = 1;
    vanz++;
    //checked FK & RK
  }    
  
  //H- -> W- + h0/H0  
  
  if(flW.ison()) {
    for(short int i=31; i<33;i++) {
      Flavour flav = Flavour(kf::code(i)); 
      
      if(flav.ison() && flHmin.ison()) {
	
	v[vanz].in[0] = flHmin;
	v[vanz].in[1] = flW;  
	v[vanz].in[2] = flav;

	//sollte so stimmen !
	
	kcpl0 = -(M_I/num_2)*g2*K_A_M(i-31,0);
	kcpl1 = kcpl0;

	v[vanz].cpl[0]  = kcpl0.Value();
	v[vanz].cpl[1]  = kcpl1.Value();
	v[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();	  	    		
	v[vanz].cpl[2]  = 0.;v[vanz].cpl[3]  = 0.;

	v[vanz].ncf   = 1;
	v[vanz].Color = new Color_Function; 
	
	v[vanz].Color->type = cf::None; 
 
	v[vanz].nlf     = 1;
	v[vanz].Lorentz = new Lorentz_Function; 
	
	v[vanz].Lorentz->type       = lf::SSV;     
	v[vanz].Lorentz->SetParticleArg(0,2,1);     	
	
	v[vanz].on      = 1;
	vanz++;
	//checked FK & RK
      }
    }
    
    //H- -> W- +  A0
    
    if(flHmin.ison() && flA0.ison()) {
      v[vanz].in[0] = flHmin;
      v[vanz].in[1] = flW;    
      v[vanz].in[2] = flA0;

      //sollte so stimmen !

      kcpl0 = -g2/num_2;
      kcpl1 = kcpl0;

      v[vanz].cpl[0]  = kcpl0.Value();
      v[vanz].cpl[1]  = kcpl1.Value();
      v[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();	  	    		
      v[vanz].cpl[2]  = 0.;v[vanz].cpl[3]  = 0.;
      
      v[vanz].ncf   = 1;
      v[vanz].Color = new Color_Function; 
      
      v[vanz].Color->type = cf::None; 
 
      v[vanz].nlf     = 1;
      v[vanz].Lorentz = new Lorentz_Function; 
      
      v[vanz].Lorentz->type       = lf::SSV;     
      v[vanz].Lorentz->SetParticleArg(0,2,1);     	
      
      v[vanz].on      = 1;
      vanz++;
      //checked FK & RK
    }
  } 
}


inline Kabbala Model_Higgs::K_yuk(Flavour fl) {
  return Kabbala(string("M_{"+fl.texname()+"}"),fl.yuk());
}

Kabbala Model_Higgs::K_yuk_sign(Flavour fl) {
  char hi[3];
  sprintf(hi,"%i",fl.get_mass_sign());
  return Kabbala(string(hi),fl.get_mass_sign());
}

inline Kabbala Model_Higgs::K_v0() {
  return Kabbala(string("v_0"),SpHiggs.v_0());
}  

inline Kabbala Model_Higgs::K_v1() {
  return Kabbala(string("v_1"),SpHiggs.v_1());
}  

inline Kabbala Model_Higgs::K_v2() {
  return Kabbala(string("v_2"),SpHiggs.v_2());
}

inline Kabbala Model_Higgs::K_Z_H(short int i,short int j) {   
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("Z^{")+string(hi)+string(hj)+string("}_H"),SpHiggs.Z_H(i,j));
}     

inline Kabbala Model_Higgs::K_Z_R(short int i,short int j) {   
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("Z^{")+string(hi)+string(hj)+string("}_R"),SpHiggs.Z_R(i,j));
}  

inline Kabbala Model_Higgs::K_A_M(short int i,short int j) {
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("A^{")+string(hi)+string(hj)+string("}_M"),
		 SpHiggs.Z_R(0,i)*SpHiggs.Z_H(0,j)-SpHiggs.Z_R(1,i)*SpHiggs.Z_H(1,j));
}  

inline Kabbala Model_Higgs::K_A_H(short int i,short int j) {
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("A^{")+string(hi)+string(hj)+string("}_H"),
		 SpHiggs.Z_H(0,i)*SpHiggs.Z_H(0,j)-SpHiggs.Z_H(1,i)*SpHiggs.Z_H(1,j));
}  

inline Kabbala Model_Higgs::K_A_P(short int i,short int j) {
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("A^{")+string(hi)+string(hj)+string("}_P"),
		 SpHiggs.Z_R(0,i)*SpHiggs.Z_H(1,j)+SpHiggs.Z_R(1,i)*SpHiggs.Z_H(0,j));
}  

inline Kabbala Model_Higgs::K_B_R(short int i) {   
  char hi[2];
  sprintf(hi,"%i",i);
  return Kabbala(string("B^{")+string(hi)+string("}_R"),
		 SpHiggs.v_1()*SpHiggs.Z_R(0,i)-SpHiggs.v_2()*SpHiggs.Z_R(1,i));
}  
























