#include "Model_sLeptons.H"
#include <stdio.h>

using namespace AMEGIC;
using namespace APHYTOOLS;
using namespace std;

Model_sLeptons::~Model_sLeptons()
{
  if (ini_higgs==1) delete moHiggs;
}

void Model_sLeptons::Init()
{
  if (moHiggs==NULL) {
    //THDM part
    moHiggs = new Model_Higgs;
    moHiggs->Init();
    ini_higgs = 1;
  }
  // Spectrum
  if (isa!=NULL) SpsN.Interface(isa);
  else
  SpsN.Init();
  if (isa!=NULL) SpsL.Interface(isa);
  else
  SpsL.Init();

  g1     = Kabbala(string("e"),sqrt(4.*M_PI*moHiggs->Aqed()));
  g2     = Kabbala(string("e/sin\\theta_W"), g1.Value()/moHiggs->SinTW());
  PL     = Kabbala(string("P_L"),1.);
  PR     = Kabbala(string("P_R"),1.);
  M_I    = Kabbala(string("i"),Complex(0.,1.));
  root2  = Kabbala(string("\\sqrt{2}"),sqrt(2.));
  K_zero = Kabbala(string("zero"),0.);
  num_2  = Kabbala(string("2"),2.);    	
  num_4  = Kabbala(string("4"),4.);    		
}

void Model_sLeptons::Init(Model_Higgs* _moHiggs,Isajet* _isa)
{
  moHiggs = _moHiggs;
  isa = _isa;
  Init();
}

void Model_sLeptons::c_SSS(Single_Vertex* v,int& vanz)
{
  Kabbala kcpl0,kcpl1;

  //sneutrino - Higgs - sneutrino
  for (short int i=81;i<84;i++) {
    Flavour flav = Flavour(kf::code(i));
    for (short int k=31;k<33;k++) {
      Flavour flh = Flavour(kf::code(k));
      if (flh.ison() && flav.ison()) {
	v[vanz].in[0] = flav;
	v[vanz].in[1] = flh;
	v[vanz].in[2] = flav;

	kcpl0 = -M_I*g2*g2/(K_cosTW()*K_cosTW()*num_4)*K_B_R(k-31);
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
	
	v[vanz].Lorentz->type       = lf::SSS;     
	
	v[vanz].on      = 1;
	//v[vanz].on      = 0;
	vanz++;
	//checked RK & FK
      }
    }
  }

  //sneutrino - Hmin - slepton
 
  Flavour flHm = Flavour(kf::Hmin);
  if (flHm.ison()) {
    for (short int i=81;i<84;i++) {
      Flavour flav1 = Flavour(kf::code(i));
      for (short int j=71;j<77;j++) {
	Flavour flav2 =Flavour(kf::code(j));
	if(flav1.ison() && flav2.ison()){
	  v[vanz].in[0] = flav1;
	  v[vanz].in[1] = flHm;
	  v[vanz].in[2] = flav2.bar();
     
	  Kabbala K_lI = Kabbala(string("\\frac{\\m M_{")+flav2.texname()+string("}}{ v_1}\\sqrt{2}"),
				 -Flavour(kf::code(2*gen_sLep(flav2)+11)).yuk()/K_v1().Value()*sqrt(2.));

	  kcpl0 = M_I*K_Z_nue(gen_sLep(flav2),i-81)*
	    (-root2*g2*g2/num_4*(K_v1()*K_Z_H(0,0)+K_v2()*K_Z_H(1,0))*
	     K_Z_L(gen_sLep(flav2),j-71)-
	     K_Z_H(0,0)*(K_yuk(flav2)*K_lI*K_Z_L(gen_sLep(flav2),j-71)-
	      (K_l_S(gen_sLep(flav2),0)*K_Z_L(3,j-71)+K_l_S(gen_sLep(flav2),1)*K_Z_L(4,j-71)+
	       K_l_S(gen_sLep(flav2),2)*K_Z_L(5,j-71)))+
	     K_Z_H(1,0)*(K_k_S(gen_sLep(flav2),0)*K_Z_L(3,j-71)+K_k_S(gen_sLep(flav2),1)*K_Z_L(4,j-71)+
			 K_k_S(gen_sLep(flav2),2)*K_Z_L(5,j-71)-K_lI*K_h()*K_Z_L(gen_sLep(flav2)+3,j-71)));
	 
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
	  
	  v[vanz].Lorentz->type       = lf::SSS;     
	  
	  v[vanz].on      = 1;
	  vanz++;
	  
	}
      }
    }
  }
  
  //slepton - A0 - slepton
  Flavour flA0 = Flavour(kf::A0);
  if (flA0.ison()) {
    for (short int i=71;i<77;i++) {
      Flavour flav1 = Flavour(kf::code(i));
      for (short int j=71;j<77;j++) {
	Flavour flav2 =Flavour(kf::code(j));
	if(flav1.ison() && flav2.ison() && i<=j){
	  v[vanz].in[0] = flav1;
	  v[vanz].in[1] = flA0;
	  v[vanz].in[2] = flav2;

	  Kabbala K_lI = Kabbala(string("\\frac{(\\m M_{")+flav1.texname()+string("})}{ v_1}\\sqrt{2}"),
				 -Flavour(kf::code(2*gen_sLep(flav1)+11)).yuk()/(K_v1()).Value()*sqrt(2.));
     
	  Kabbala inv_root2  = Kabbala(string("\\frac{1}{\\sqrt{2}}"),1./sqrt(2.));

	  kcpl0 = (K_l_S(gen_sLep(flav1),gen_sLep(flav2))*
		    (K_Z_L(gen_sLep(flav1),j-71)*K_Z_L(gen_sLep(flav2)+3,i-71)-
		     K_Z_L(gen_sLep(flav1),i-71)*K_Z_L(gen_sLep(flav2)+3,j-71))*K_Z_H(0,0)+
		    K_k_S(gen_sLep(flav1),gen_sLep(flav2))*
		    (K_Z_L(gen_sLep(flav1),j-71)*K_Z_L(gen_sLep(flav2)+3,i-71)-
		     K_Z_L(gen_sLep(flav1),i-71)*K_Z_L(gen_sLep(flav2)+3,j-71))*K_Z_H(1,0)+
		    K_lI*K_h()*(K_Z_L(gen_sLep(flav1),i-71)*K_Z_L(gen_sLep(flav1)+3,j-71)-
				K_Z_L(gen_sLep(flav1),j-71)*K_Z_L(gen_sLep(flav1)+3,i-71))*
		    K_Z_H(1,0))*(-inv_root2);
	  
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
	  
	  v[vanz].Lorentz->type       = lf::SSS;     
	  	  
	  v[vanz].on      = 1;
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
	if(flav1.ison() && flav2.ison() && i<=j){
	  v[vanz].in[0] = flav1;
	  v[vanz].in[1] = flH;
	  v[vanz].in[2] = flav2;

	  Kabbala help = K_zero;

	  Kabbala K_lI = Kabbala(string("\\frac{\\m M_{")+flav1.texname()+string("}}{ v_1}*\\sqrt{2}"),
				 -Flavour(kf::code(2*gen_sLep(flav1)+11)).yuk()/(K_v1()).Value()*sqrt(2.));
	 
	  Kabbala fac = Kabbala(string("\\frac{1-4sin^2\\theta_W}{2sin^2\\theta_W}"),
				(1.-4.*(K_sinTW()).Value()*(K_sinTW()).Value())/
				(2.*(K_sinTW()).Value()*(K_sinTW()).Value()));

	  if (i==j) {help = Kabbala(string("1"),1.);}

	  kcpl0 = M_I*(g1*g1/(K_cosTW()*K_cosTW()*num_2)*K_B_R(k-31)*
		       (help+fac*K_Z_L(gen_sLep(flav1),i-71)*K_Z_L(gen_sLep(flav1),j-71))-
		       K_lI*K_lI*K_v1()*K_Z_R(0,k-31)*
		       (K_Z_L(gen_sLep(flav1),i-71)*K_Z_L(gen_sLep(flav1),j-71)+
			K_Z_L(gen_sLep(flav1)+3,i-71)*K_Z_L(gen_sLep(flav1)+3,j-71))-
		       K_Z_R(0,k-31)/root2*K_l_S(gen_sLep(flav1),gen_sLep(flav2))*
						  (K_Z_L(gen_sLep(flav1),j-71)*
						   K_Z_L(gen_sLep(flav2)+3,i-71)+
						   K_Z_L(gen_sLep(flav1),i-71)*
						   K_Z_L(gen_sLep(flav2)+3,j-71))+
		       K_Z_R(1,k-31)/root2*K_k_S(gen_sLep(flav1),gen_sLep(flav2))*
						  (K_Z_L(gen_sLep(flav1),j-71)*
						   K_Z_L(gen_sLep(flav2)+3,i-71)+
						   K_Z_L(gen_sLep(flav1),i-71)*
						   K_Z_L(gen_sLep(flav2)+3,j-71))-
		       K_lI*K_Z_R(1,k-31)*K_h()/root2*(K_Z_L(gen_sLep(flav1),i-71)*
						       K_Z_L(gen_sLep(flav1)+3,j-71)+
						       K_Z_L(gen_sLep(flav1),j-71)*
						       K_Z_L(gen_sLep(flav1)+3,i-71)));
	  
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
	  
	  v[vanz].Lorentz->type       = lf::SSS;     
	  
	  v[vanz].on      = 1;
	  //v[vanz].on      = 0;
	  vanz++;
  	}
      }
    }
  }
}

void Model_sLeptons::c_SSV(Single_Vertex* v,int& vanz)
{
  Kabbala kcpl0,kcpl1;

  //slepton - Photon - slepton
  Flavour flPh = Flavour(kf::photon);
  if (flPh.ison()) {
    for (short int i=71;i<77;i++) {
      Flavour flav = Flavour(kf::code(i));
      Kabbala charge1 = Kabbala(string("Q_{")+flav.texname()+string("}"),flav.charge());
      if (flav.ison()) {
	v[vanz].in[0] = flav;
	v[vanz].in[1] = flPh;
	v[vanz].in[2] = flav;
	
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
      
	v[vanz].Lorentz->type       = lf::SSV;     
	v[vanz].Lorentz->SetParticleArg(0,2,1);     	

	v[vanz].on      = 1;
	vanz++;
      }
    }
  }
 
 //slepton - Z - slepton
  Flavour flZ = Flavour(kf::Z);
  if (flZ.ison()) {
    for (short int i=71;i<77;i++) {
      Flavour flav1 = Flavour(kf::code(i));
      for (short int j=i;j<77;j++) {
	Flavour flav2 = Flavour(kf::code(j));
	if (flav1.ison() && flav2.ison()) {
	  
	  Kabbala help = K_zero;
	  if(i==j) help = K_sinTW()*K_sinTW();
	  
	  v[vanz].in[0] = flav1;
	  v[vanz].in[1] = flZ;
	  v[vanz].in[2] = flav2;
	  	  
	  kcpl0 = M_I*g2/K_cosTW()*
	    ((K_Z_L(0,j-71)*K_Z_L(0,i-71)+
	      K_Z_L(1,j-71)*K_Z_L(1,i-71)+
	      K_Z_L(2,j-71)*K_Z_L(2,i-71))/num_2-help);
	  	
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
	  
	  v[vanz].Lorentz->type       = lf::SSV;     
	  v[vanz].Lorentz->SetParticleArg(0,2,1);     	
	  
	  v[vanz].on      = 1;
	  vanz++;
	  //checked SS
	}
      }

    }

    //sneutrino - Z - sneutrino
    for (short int i=81;i<84;i++) {
      Flavour flav = Flavour(kf::code(i));
      if (flav.ison()) {
 	  v[vanz].in[0] = flav;
	  v[vanz].in[1] = flZ;
	  v[vanz].in[2] = flav;

	  //changed sign
	  kcpl0 = -M_I*g2/(K_cosTW()*num_2);
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
	  
	  v[vanz].Lorentz->type       = lf::SSV;     
	  v[vanz].Lorentz->SetParticleArg(0,2,1);     	
	  
	  v[vanz].on      = 1;
	  vanz++;
	  //checked SS
      }
    }
    
  }
  
  //check for summing convention !!!
  //sneutrino - W - slepton
  Flavour flW = Flavour(kf::W);
  if (flW.ison()) {
      for (short int i=81;i<84;i++) {
      Flavour flav1 = Flavour(kf::code(i));
      for (short int j=71;j<77;j++) {
	Flavour flav2 =Flavour(kf::code(j));
	if(flav1.ison() && flav2.ison()){
	  v[vanz].in[0] = flav1;
	  v[vanz].in[1] = flW;
	  v[vanz].in[2] = flav2.bar();
	
	  kcpl0 = -M_I*g2/root2*K_Z_nue(gen_sLep(flav2),i-81)*
	    K_Z_L(gen_sLep(flav2),j-71);
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
	  
	  v[vanz].Lorentz->type       = lf::SSV;     
	  v[vanz].Lorentz->SetParticleArg(0,2,1);     	
	  
	  v[vanz].on      = 1;
	  vanz++;
  	}
      }
    }
  }
}

inline Kabbala Model_sLeptons::K_Z_nue(short int i,short int j)       
{   
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("Z^{")+string(hi)+string(hj)+string("}_\\nu"),SpsN.Znue(i,j));
}  
inline Kabbala Model_sLeptons::K_Z_L(short int i,short int j)       
{   
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("Z^{")+string(hi)+string(hj)+string("}_L"),SpsL.Z_L(i,j));
}   
inline Kabbala Model_sLeptons::K_l_S(short int i,short int j)
{   
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("l^{")+string(hi)+string(hj)+string("}_S"),SpsL.l_S(i,j));
}  
inline Kabbala Model_sLeptons::K_k_S(short int i,short int j)
{   
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("k^{")+string(hi)+string(hj)+string("}_S"),SpsL.k_S(i,j));
}  
inline int Model_sLeptons::gen_sLep(Flavour fl)
{
  int gen_sL;

  if (fl.kfcode() == 71 || fl.kfcode() == 74)
    gen_sL = 0;
  if (fl.kfcode() == 72 || fl.kfcode() == 75)
    gen_sL = 1;
  if (fl.kfcode() == 73 || fl.kfcode() == 76)
    gen_sL = 2;

  return gen_sL;
}
