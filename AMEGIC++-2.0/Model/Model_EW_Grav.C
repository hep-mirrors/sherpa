#include <stdio.h>
#include "Model_EW_Grav.H"
#include "Running_AlphaQED.H"
#include "Run_Parameter.H"
#include "MathTools.H"
#include "Vector.H"

using namespace AMEGIC;
using namespace APHYTOOLS;
using namespace AMATOOLS;
using namespace AORGTOOLS;
using namespace std;

void Model_EW_Grav::Init()
{
  // Initialize Yukawa masses
  SpEW.FillYukawas();
  //Couplings
  CplEW.Init();
  CplLED.Init();

  g1    = Kabbala(string("g_1"),sqrt(4.*M_PI*Aqed()));
  g2    = Kabbala(string("g_1/sin\\theta_W"), g1.Value()/SinTW());
  PL    = Kabbala(string("P_L"),1.);
  PR    = Kabbala(string("P_R"),1.);
  M_I   = Kabbala(string("i"),Complex(0.,1.));
  root2 = Kabbala(string("\\sqrt{2}"),sqrt(2.));
  //vev   = Kabbala(string("v_{EW}"),CplEW.VEV());
  vev   = Kabbala(string("v_{EW}"),2.*SinTW()*CosTW()*Flavour(kf::Z).Mass()/g1.Value());
  cout<<"VeVs: "<<CplEW.VEV()<<" =?= "<<vev.Value()<<endl;
  kap   = Kabbala(string("kappa"),CplLED.Kappa());
  om    = Kabbala(string("omega"),CplLED.Omega());
  num2  = Kabbala(string("2"),2.);
  num4  = Kabbala(string("4"),4.);
  num15 = Kabbala(string("1.5"),1.5);

}

void Model_EW_Grav::c_FFT(Single_Vertex* v,int& vanz)
{
  Flavour flgraviton(kf::graviton);
  Flavour flgs(kf::gscalar);
  for (short int i=1;i<17;i++) {
    if (i==7) i=11;
    Flavour flav1 = Flavour(kf::code(i));
	
    if (flav1.IsOn()) {
      for (short int j=i;j<17;j++) {
	if (j==7) j=11;	
	Flavour flav2 = Flavour(kf::code(j));
	
	if (flav2.IsOn()) {
	  if (flav1==flav2) {
	    
	    if (flgraviton.IsOn()) {
	      
	      Kabbala mf = Kabbala(string("M_{")+flav1.TexName()+string("}"),flav1.Mass());
	      Kabbala kcpl0,kcpl1;
	      kcpl0 = -M_I*kap/num4;
	      kcpl1 = num2*mf;
	      //kcpl0=Kabbala("",5520);
	      //kcpl1=Kabbala("",5521);
	      
	      v[vanz].in[0] = flav1;
	      v[vanz].in[1] = Flavour(kf::graviton);
	      v[vanz].in[2] = flav2;
	      v[vanz].cpl[0]  = kcpl0.Value();
	      v[vanz].cpl[1]  = v[vanz].cpl[0];
	      v[vanz].Str     = (kcpl0*PR+kcpl0*PL).String();
	      v[vanz].cpl[2]  = kcpl1.Value();
	      v[vanz].cpl[3]  = 0.;
	      
	      v[vanz].ncf   = 1;
	      v[vanz].Color = new Color_Function; 
	      
	      if (flav1.Strong()) {
		v[vanz].Color->type       = cf::D;     
		v[vanz].Color->SetParticleArg(0,2);     
		v[vanz].Color->SetStringArg('0','2');     
	      }
	      else v[vanz].Color->type = cf::None; 
	      
	      v[vanz].nlf     = 1;
	      v[vanz].Lorentz = new Lorentz_Function; 
	      
	      v[vanz].Lorentz->type       = lf::FFT;     
	      v[vanz].Lorentz->SetParticleArg(1);     

	      v[vanz].on      = 1;
	      vanz++;	      
	    }
	    //scalar graviton mode
	    if (flgs.IsOn()) {
	      
	      Kabbala mf = Kabbala(string("M_{")+flav1.TexName()+string("}"),flav1.Mass());
	      Kabbala kcpl0,kcpl1;
	      kcpl0 = M_I*om*kap*Kabbala(string("3/4"),.75);
	      kcpl1 = mf*Kabbala(string("8/3"),8./3.);
	      //kcpl0=Kabbala("",5500);
	      //kcpl1=Kabbala("",5501);
	     
	      v[vanz].in[0] = flav1;
	      v[vanz].in[1] = Flavour(kf::gscalar);
	      v[vanz].in[2] = flav2;
	      v[vanz].cpl[0] = kcpl0.Value();
	      v[vanz].cpl[1] = kcpl0.Value();
	      v[vanz].Str     = (kcpl0*PR+kcpl0*PL).String();
	      v[vanz].cpl[2] = kcpl1.Value(); 
	      v[vanz].cpl[3]  = 0.;
	
	      v[vanz].ncf   = 1;
	      v[vanz].Color = new Color_Function; 
	      
	      if (flav1.Strong()) {
		v[vanz].Color->type       = cf::D;     
		v[vanz].Color->SetParticleArg(0,2);     
		v[vanz].Color->SetStringArg('0','2');     
	      }
	      else v[vanz].Color->type = cf::None; 
	      
	      v[vanz].nlf     = 1;
	      v[vanz].Lorentz = new Lorentz_Function; 
	      
	      v[vanz].Lorentz->type       = lf::FFGS;     
	      //cout<< v[vanz].in[0]<<v[vanz].in[1]<<v[vanz].in[2]<<endl;
	      v[vanz].on     = 1;
	      vanz++;
	    }
	  }   
     	}
      }
    }
  }
}


void Model_EW_Grav::c_FFVT(Single_Vertex* v,int& vanz)
{
  //return;
  Flavour flphoton(kf::photon);
  Flavour flZ(kf::Z);
  Flavour flW(kf::W);
  Flavour flgraviton(kf::graviton);
  Flavour flgs(kf::gscalar);
  for (short int i=1;i<17;i++) {
    if (i==7) i=11;
    Flavour flav1 = Flavour(kf::code(i));
    Kabbala charge1 = Kabbala(string("Q_{")+flav1.TexName()+string("}"),flav1.Charge());
    Kabbala isoweak1 = Kabbala(string("T_{")+flav1.TexName()+string("}"),flav1.IsoWeak());
    
    if (flav1.IsOn()) {
      for (short int j=i;j<17;j++) {
	if (j==7) j=11;	
	Flavour flav2 = Flavour(kf::code(j));
	Kabbala charge2 = Kabbala(string("Q_{")+flav2.TexName()+string("}"),flav2.Charge());
	Kabbala isoweak2 = Kabbala(string("T_{")+ flav2.TexName()+string("}"),flav2.IsoWeak());	
	      
	if (flav2.IsOn()) {
	  if (flav1==flav2) {
	    //photon + graviton
	    if (flphoton.IsOn()&&flgraviton.IsOn()) {
	      
	      Kabbala kcpl0,kcpl1;
	      kcpl0 = g1*M_I*charge1*kap/num2;
	      kcpl1 = kcpl0;
	      //kcpl0=Kabbala("",55120);
	      //kcpl1=Kabbala("",55121);
	      
	      if (!AMATOOLS::IsZero(charge1.Value())) {
		v[vanz].nleg     = 4;
		v[vanz].in[0] = flav1;
		v[vanz].in[1] = Flavour(kf::photon);
		v[vanz].in[2] = flav2;
		v[vanz].in[3] = flgraviton;
		v[vanz].cpl[0]  = kcpl0.Value();
		v[vanz].cpl[1]  = v[vanz].cpl[0];
		v[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
		v[vanz].cpl[2]  = 0.;v[vanz].cpl[3]  = 0.;
	
		v[vanz].ncf   = 1;
		v[vanz].Color = new Color_Function; 
		
		if (flav1.Strong()) {
		  v[vanz].Color->type       = cf::D;     
		  v[vanz].Color->SetParticleArg(0,2);     
		  v[vanz].Color->SetStringArg('0','2');     
		}
		else v[vanz].Color->type = cf::None; 
		
		v[vanz].nlf     = 1;
		v[vanz].Lorentz = new Lorentz_Function; 
		
		v[vanz].Lorentz->type       = lf::FFVT;     
		v[vanz].Lorentz->SetParticleArg(1,3);     

		v[vanz].on      = 1;
		vanz++;
	      }
	    }
	    //Z +graviton
	    if (flZ.IsOn()&&flgraviton.IsOn()) {
	      
	      Kabbala kcpl0,kcpl1;
	      kcpl0 = -M_I/K_cosTW()*charge1*K_sinTW()*K_sinTW()*g2*kap/num2;
	      kcpl1 = M_I/K_cosTW()*(isoweak1-charge1*K_sinTW()*K_sinTW())*g2*kap/num2;
	      //kcpl0=Kabbala("",55120);
	      //kcpl1=Kabbala("",55121);

	      v[vanz].nleg     = 4;
	      v[vanz].in[0] = flav1;
	      v[vanz].in[1] = Flavour(kf::Z);
	      v[vanz].in[2] = flav2;
	      v[vanz].in[3] = flgraviton;
	      v[vanz].cpl[0] = kcpl0.Value();
	      v[vanz].cpl[1] = kcpl1.Value();
	      v[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	      v[vanz].cpl[2] = 0.;v[vanz].cpl[3]  = 0.;
	
	      v[vanz].ncf   = 1;
	      v[vanz].Color = new Color_Function; 
	      
	      if (flav1.Strong()) {
		v[vanz].Color->type       = cf::D;     
		v[vanz].Color->SetParticleArg(0,2);     
		v[vanz].Color->SetStringArg('0','2');     
	      }
	      else v[vanz].Color->type = cf::None; 
	      
	      v[vanz].nlf     = 1;
	      v[vanz].Lorentz = new Lorentz_Function; 
	      
	      v[vanz].Lorentz->type       = lf::FFVT;     
	      v[vanz].Lorentz->SetParticleArg(1,3);     

	      v[vanz].on     = 1;
	      vanz++;
	    }
	    //photon + gscalar
	    if (flphoton.IsOn()&&flgs.IsOn()) {
	      
	      Kabbala kcpl0,kcpl1;
	      kcpl0 = -g1*M_I*charge1*om*kap*num15;
	      kcpl1 = kcpl0;
	      //kcpl0=Kabbala("",55100);
	      //kcpl1=Kabbala("",55101);
	      
	      if (!AMATOOLS::IsZero(charge1.Value())) {
		v[vanz].nleg     = 4;
		v[vanz].in[0] = flav1;
		v[vanz].in[1] = Flavour(kf::photon);
		v[vanz].in[2] = flav2;
		v[vanz].in[3] = flgs;
		v[vanz].cpl[0]  = kcpl0.Value();
		v[vanz].cpl[1]  = v[vanz].cpl[0];
		v[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
		v[vanz].cpl[2]  = 0.;v[vanz].cpl[3]  = 0.;
	
		v[vanz].ncf   = 1;
		v[vanz].Color = new Color_Function; 
		
		if (flav1.Strong()) {
		  v[vanz].Color->type       = cf::D;     
		  v[vanz].Color->SetParticleArg(0,2);     
		  v[vanz].Color->SetStringArg('0','2');     
		}
		else v[vanz].Color->type = cf::None; 
		
		v[vanz].nlf     = 1;
		v[vanz].Lorentz = new Lorentz_Function; 
		
		v[vanz].Lorentz->type       = lf::FFVGS;     
		v[vanz].Lorentz->SetParticleArg(1);     

		v[vanz].on      = 1;
		vanz++;
	      }
	    }
	    //Z +gscalar
	    if (flZ.IsOn()&&flgs.IsOn()) {
	      
	      Kabbala kcpl0,kcpl1;
	      kcpl0 = M_I/K_cosTW()*charge1*K_sinTW()*K_sinTW()*g2*om*kap*num15;
	      kcpl1 = -M_I/K_cosTW()*(isoweak1-charge1*K_sinTW()*K_sinTW())*g2*om*kap*num15;
	      //kcpl0=Kabbala("",55100);
	      //kcpl1=Kabbala("",55101);

	      v[vanz].nleg     = 4;
	      v[vanz].in[0] = flav1;
	      v[vanz].in[1] = Flavour(kf::Z);
	      v[vanz].in[2] = flav2;
	      v[vanz].in[3] = flgs;
	      v[vanz].cpl[0] = kcpl0.Value();
	      v[vanz].cpl[1] = kcpl1.Value();
	      v[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	      v[vanz].cpl[2] = 0.;v[vanz].cpl[3]  = 0.;
	
	      v[vanz].ncf   = 1;
	      v[vanz].Color = new Color_Function; 
	      
	      if (flav1.Strong()) {
		v[vanz].Color->type       = cf::D;     
		v[vanz].Color->SetParticleArg(0,2);     
		v[vanz].Color->SetStringArg('0','2');     
	      }
	      else v[vanz].Color->type = cf::None; 
	      
	      v[vanz].nlf     = 1;
	      v[vanz].Lorentz = new Lorentz_Function; 
	      
	      v[vanz].Lorentz->type       = lf::FFVGS;     
	      v[vanz].Lorentz->SetParticleArg(1);     

	      v[vanz].on     = 1;
	      vanz++;
	    }
	  }
	  //W + graviton
	  if (flW.IsOn()&&flgraviton.IsOn()) {
	    short int hit = 1;
	    Kabbala kcpl0,kcpl1;
	    kcpl0 = Kabbala(string("zero"),0.);
	    kcpl1 = Kabbala(string("1"),0.);

	    if (!((flav1.IsDowntype() && flav2.IsUptype()) ||
                  (flav2.IsDowntype() && flav1.IsUptype()))) hit = 0;
	    if ((flav1.IsLepton() && !flav2.IsLepton()) ||
		(flav1.IsQuark() && !flav2.IsQuark()) ) hit = 0;
	    if (hit==1) {
	      if (flav1.IsDowntype() && i>10 && j==i+1) 
		kcpl1 = M_I/root2*g2*kap/num2;
	      if (i<7 && j<7) {
		if (flav1.IsDowntype())
		kcpl1 = M_I/root2*g2*K_CKM((i-1)/2,j/2-1)*kap/num2;
		else 	    
		kcpl1 = M_I/root2*g2*K_CKM(i/2-1,(j-1)/2)*kap/num2;		
	      }
	      //cout<<"W-Decay: "<<flav1<<";"<<flav2<<"  :  "<<c1<<endl;
	      if (!AMATOOLS::IsZero(kcpl1.Value()/kap.Value())) {
		v[vanz].nleg     = 4;
		v[vanz].in[1] = Flavour(kf::W);
		if (flav1.IsDowntype()) {
		  v[vanz].in[0] = flav1;
		  v[vanz].in[2] = flav2;
		}
		else {
		  v[vanz].in[0] = flav2;
		  v[vanz].in[2] = flav1;
		}
		v[vanz].in[3] = flgraviton;
		
		v[vanz].cpl[0]  = kcpl0.Value();
		v[vanz].cpl[1]  = kcpl1.Value();
		v[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
		v[vanz].cpl[2]  = 0.;v[vanz].cpl[3]  = 0.;
	
		v[vanz].ncf   = 1;
		v[vanz].Color = new Color_Function; 
		
		if (flav1.Strong()) {
		  v[vanz].Color->type       = cf::D;     
		  v[vanz].Color->SetParticleArg(0,2);     
		  v[vanz].Color->SetStringArg('0','2');     
		}
		else v[vanz].Color->type = cf::None; 
		
		v[vanz].nlf     = 1;
		v[vanz].Lorentz = new Lorentz_Function; 
		
		v[vanz].Lorentz->type       = lf::FFVT;     
		v[vanz].Lorentz->SetParticleArg(1,3);     

		v[vanz].on      = 1;
		vanz++;
	      }
	    }
	  }
	  //W + gscalar
	  if (flW.IsOn()&&flgs.IsOn()) {
	    short int hit = 1;
	    Kabbala kcpl0,kcpl1;
	    kcpl0 = Kabbala(string("zero"),0.);
	    kcpl1 = Kabbala(string("1"),0.);

	    if (!((flav1.IsDowntype() && flav2.IsUptype()) ||
                  (flav2.IsDowntype() && flav1.IsUptype()))) hit = 0;
	    if ((flav1.IsLepton() && !flav2.IsLepton()) ||
		(flav1.IsQuark() && !flav2.IsQuark()) ) hit = 0;
	    if (hit==1) {
	      if (flav1.IsDowntype() && i>10 && j==i+1) 
		kcpl1 = -M_I/root2*g2*om*kap*num15;
	      if (i<7 && j<7) {
		if (flav1.IsDowntype())
		kcpl1 = -M_I/root2*g2*K_CKM((i-1)/2,j/2-1)*om*kap*num15;
		else 	    
		kcpl1 = -M_I/root2*g2*K_CKM(i/2-1,(j-1)/2)*om*kap*num15;		
	      }
	      //cout<<"W-Decay: "<<flav1<<";"<<flav2<<"  :  "<<c1<<endl;
	      if (!AMATOOLS::IsZero(kcpl1.Value()/kap.Value())) {
		v[vanz].nleg     = 4;
		v[vanz].in[1] = Flavour(kf::W);
		if (flav1.IsDowntype()) {
		  v[vanz].in[0] = flav1;
		  v[vanz].in[2] = flav2;
		}
		else {
		  v[vanz].in[0] = flav2;
		  v[vanz].in[2] = flav1;
		}
		v[vanz].in[3] = flgs;
		
		v[vanz].cpl[0]  = kcpl0.Value();
		v[vanz].cpl[1]  = kcpl1.Value();
		v[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
		v[vanz].cpl[2]  = 0.;v[vanz].cpl[3]  = 0.;
	
		v[vanz].ncf   = 1;
		v[vanz].Color = new Color_Function; 
		
		if (flav1.Strong()) {
		  v[vanz].Color->type       = cf::D;     
		  v[vanz].Color->SetParticleArg(0,2);     
		  v[vanz].Color->SetStringArg('0','2');     
		}
		else v[vanz].Color->type = cf::None; 
		
		v[vanz].nlf     = 1;
		v[vanz].Lorentz = new Lorentz_Function; 
		
		v[vanz].Lorentz->type       = lf::FFVGS;     
		v[vanz].Lorentz->SetParticleArg(1);     

		v[vanz].on      = 1;
		vanz++;
	      }
	    }
	  }
	}
      }
    }
  }
}

void Model_EW_Grav::c_VVT(Single_Vertex* v,int& vanz)
{
  Flavour flgraviton(kf::graviton);
  Flavour flgs(kf::gscalar);
  Kabbala kcpl0,kcpl1;  
  Flavour flav(kf::W);

  if (flgraviton.IsOn()){  
    // W graviton W
    if (flav.IsOn()) {
      v[vanz].in[0] = flav;
      v[vanz].in[1] = flgraviton;
      v[vanz].in[2] = flav;
      
      kcpl0 = -M_I*kap;
      kcpl1 = Kabbala(string("\\sqr(M_{")+flav.TexName()+string("})"),sqr(flav.Mass()));
      //kcpl0=Kabbala("",1120);
      //kcpl1=Kabbala("",1121);
      
      v[vanz].cpl[0]  = kcpl0.Value();
      v[vanz].cpl[1]  = kcpl1.Value();
      v[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
      v[vanz].cpl[2]  = 0.;v[vanz].cpl[3]  = 0.;

      v[vanz].ncf   = 1;
      v[vanz].Color = new Color_Function; 
      
      v[vanz].Color->type       = cf::None;     
      
      v[vanz].nlf     = 1;
      v[vanz].Lorentz = new Lorentz_Function; 
    
      v[vanz].Lorentz->type = lf::VVT;     
      v[vanz].Lorentz->SetParticleArg(0,2,1);     

      v[vanz].on      = 1;
      vanz++;
    }

    flav = Flavour(kf::Z);
    // Z graviton Z
    if (flav.IsOn()) {
      v[vanz].in[0] = flav;
      v[vanz].in[1] = flgraviton;
      v[vanz].in[2] = flav;
      
      kcpl0 = -M_I*kap;
      kcpl1 = Kabbala(string("\\sqr(M_{")+flav.TexName()+string("})"),sqr(flav.Mass()));
    //kcpl0=Kabbala("",1120);
    //kcpl1=Kabbala("",1121);
      
      v[vanz].cpl[0]  = kcpl0.Value();
      v[vanz].cpl[1]  = kcpl1.Value();
      v[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
      v[vanz].cpl[2]  = 0.;v[vanz].cpl[3]  = 0.;

      v[vanz].ncf   = 1;
      v[vanz].Color = new Color_Function; 
      
      v[vanz].Color->type       = cf::None;     
      
      v[vanz].nlf     = 1;
      v[vanz].Lorentz = new Lorentz_Function; 
    
      v[vanz].Lorentz->type = lf::VVT;     
      v[vanz].Lorentz->SetParticleArg(0,2,1);     

      v[vanz].on      = 1;
      vanz++;
    }
    flav = Flavour(kf::photon);
    // photon graviton photon
    if (flav.IsOn()) {
      v[vanz].in[0] = flav;
      v[vanz].in[1] = flgraviton;
      v[vanz].in[2] = flav;
    
      kcpl0 = -M_I*kap;
      kcpl1 = Kabbala("",0);
      //kcpl0=Kabbala("",1120);
      //kcpl1=Kabbala("",1121);

      v[vanz].cpl[0]  = kcpl0.Value();
      v[vanz].cpl[1]  = kcpl1.Value();
      v[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
      v[vanz].cpl[2]  = 0.;v[vanz].cpl[3]  = 0.;
      
      v[vanz].ncf   = 1;
      v[vanz].Color = new Color_Function; 
    
      v[vanz].Color->type       = cf::None;     

      v[vanz].nlf     = 1;
      v[vanz].Lorentz = new Lorentz_Function; 

      v[vanz].Lorentz->type = lf::VVT;     
      v[vanz].Lorentz->SetParticleArg(0,2,1);     
      
      v[vanz].on      = 1;
      vanz++;
    }
  }
  if (!flgs.IsOn()) return;
  
  flav = Flavour(kf::W);
  // W gscalar W
  if (flav.IsOn()) {
    v[vanz].in[0] = flav;
    v[vanz].in[1] = flgs;
    v[vanz].in[2] = flav;
    
    Kabbala ma = Kabbala(string("\\sqr(M_{")+flav.TexName()+string("})"),sqr(flav.Mass()));
    kcpl0 = M_I*om*kap;
    kcpl1 = ma;
    
    v[vanz].cpl[0]  = kcpl0.Value();
    v[vanz].cpl[1]  = kcpl1.Value();
    v[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
    v[vanz].cpl[2]  = 0.;v[vanz].cpl[3]  = 0.;

    v[vanz].ncf   = 1;
    v[vanz].Color = new Color_Function; 
    
    v[vanz].Color->type       = cf::None;     

    v[vanz].nlf     = 1;
    v[vanz].Lorentz = new Lorentz_Function; 

    v[vanz].Lorentz->type = lf::VVGS;     
    v[vanz].Lorentz->SetParticleArg(0,2,1);     

    v[vanz].on      = 1;
    vanz++;
  }

  flav = Flavour(kf::Z);
  // Z gscalar Z
  if (flav.IsOn()) {
    v[vanz].in[0] = flav;
    v[vanz].in[1] = flgs;
    v[vanz].in[2] = flav;
    
    Kabbala ma = Kabbala(string("\\sqr(M_{")+flav.TexName()+string("})"),sqr(flav.Mass()));
    kcpl0 = M_I*om*kap;
    kcpl1 = ma;
    //kcpl0=Kabbala("",1100);
    //kcpl1=Kabbala("",1101);

    v[vanz].cpl[0]  = kcpl0.Value();
    v[vanz].cpl[1]  = kcpl1.Value();
    v[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
    v[vanz].cpl[2]  = 0.;v[vanz].cpl[3]  = 0.;

    v[vanz].ncf   = 1;
    v[vanz].Color = new Color_Function; 
    
    v[vanz].Color->type       = cf::None;     

    v[vanz].nlf     = 1;
    v[vanz].Lorentz = new Lorentz_Function; 

    v[vanz].Lorentz->type = lf::VVGS;     
    v[vanz].Lorentz->SetParticleArg(0,2,1);     

    v[vanz].on      = 1;
    vanz++;
  }
  flav = Flavour(kf::photon);
  // photon gscalar photon
  if (flav.IsOn()) {
    v[vanz].in[0] = flav;
    v[vanz].in[1] = flgs;
    v[vanz].in[2] = flav;
    
    kcpl0 = M_I*om*kap;
    kcpl1 = Kabbala();
    //kcpl0=Kabbala("",1100);
    //kcpl1=Kabbala("",1101);

    v[vanz].cpl[0]  = kcpl0.Value();
    v[vanz].cpl[1]  = kcpl1.Value();
    v[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
    v[vanz].cpl[2]  = 0.;v[vanz].cpl[3]  = 0.;

    v[vanz].ncf   = 1;
    v[vanz].Color = new Color_Function; 
    
    v[vanz].Color->type       = cf::None;     

    v[vanz].nlf     = 1;
    v[vanz].Lorentz = new Lorentz_Function; 

    v[vanz].Lorentz->type = lf::VVGS;     
    v[vanz].Lorentz->SetParticleArg(0,2,1);     

    v[vanz].on      = 1;
    vanz++;
  }
}

void Model_EW_Grav::c_VVVT(Single_Vertex* v,int& vanz)
{
  Flavour flgraviton(kf::graviton);
  Flavour flav(kf::W);
  Kabbala kcpl0,kcpl1,kcpl0_1,kcpl1_1,charge;
  charge = Kabbala(string("Q_{")+flav.TexName()+string("}"),flav.Charge());

  if (!flav.IsOn()) return;
  if (flgraviton.IsOn()){
    // photon WW graviton
    v[vanz].nleg     = 4;
    v[vanz].in[0] = flav;
    v[vanz].in[1] = Flavour(kf::photon);
    v[vanz].in[2] = flav;
    v[vanz].in[3] = flgraviton;
   
    kcpl0 = -M_I*g1*charge*kap;
    kcpl1 = kcpl0;
    //kcpl0=Kabbala("",11120);
    //kcpl1=Kabbala("",11121);

    v[vanz].cpl[0]  = kcpl0.Value();
    v[vanz].cpl[1]  = v[vanz].cpl[0];
    v[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
    v[vanz].cpl[2]  = 0.;v[vanz].cpl[3]  = 0.;
    
    v[vanz].ncf   = 1;
    v[vanz].Color = new Color_Function; 
    
    v[vanz].Color->type       = cf::None;     
    
    v[vanz].nlf     = 1;
    v[vanz].Lorentz = new Lorentz_Function; 

    v[vanz].Lorentz->type = lf::VVVT;     
    v[vanz].Lorentz->SetParticleArg(0,1,2,3);     
    
    v[vanz].on      = 1;
    vanz++;
    
    // ZWW graviton
    v[vanz].nleg     = 4;
    v[vanz].in[0] = flav;
    v[vanz].in[1] = Flavour(kf::Z);
    v[vanz].in[2] = flav;
    v[vanz].in[3] = flgraviton;
    
    kcpl0 = -M_I*g2*charge*K_cosTW()*kap;
    kcpl1 = kcpl0;
    //kcpl0=Kabbala("",11120);
    //kcpl1=Kabbala("",11121);
  
    v[vanz].cpl[0]  = kcpl0.Value();
    v[vanz].cpl[1]  = kcpl1.Value();
    v[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
    v[vanz].cpl[2]  = 0.;v[vanz].cpl[3]  = 0.;
    
    v[vanz].ncf   = 1;
    v[vanz].Color = new Color_Function; 
    
    v[vanz].Color->type       = cf::None;     

    v[vanz].nlf     = 1;
    v[vanz].Lorentz = new Lorentz_Function; 
    
    v[vanz].Lorentz->type = lf::VVVT;     
    v[vanz].Lorentz->SetParticleArg(0,1,2,3);     
    
    v[vanz].on      = 1;
    vanz++;
  }
}


void Model_EW_Grav::c_SST(Single_Vertex* v,int& vanz)
{
  Kabbala kcpl0,kcpl1;

  Flavour flgraviton(kf::graviton);
  Flavour flgs(kf::gscalar);
  Flavour flh = Flavour(kf::h);

 if (flh.IsOn()&&flgraviton.IsOn()) {  
    v[vanz].in[0] = flh;
    v[vanz].in[1] = flgraviton;
    v[vanz].in[2] = flh;

    kcpl0 = -M_I*kap;
    kcpl1 = Kabbala(string("\\sqr(M_{")+flh.TexName()+string("})"),sqr(flh.Yuk()));
    //kcpl0=Kabbala("",320);
    //kcpl1=Kabbala("",321);
    
    v[vanz].cpl[0]  = kcpl0.Value();
    v[vanz].cpl[1]  = kcpl1.Value();
    v[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
    v[vanz].cpl[2]  = 0.;v[vanz].cpl[3]  = 0.;

    v[vanz].ncf   = 1;
    v[vanz].Color = new Color_Function; 
    
    v[vanz].Color->type   = cf::None;     
    
    v[vanz].nlf     = 1;
    v[vanz].Lorentz = new Lorentz_Function; 
    
    v[vanz].Lorentz->type = lf::SST;     
    v[vanz].Lorentz->SetParticleArg(0,2,1);     	

    v[vanz].on      = 1;
    vanz++;
  }

 if (flh.IsOn()&&flgs.IsOn()) {  
    v[vanz].in[0] = flh;
    v[vanz].in[1] = flgs;
    v[vanz].in[2] = flh;

    kcpl0 = -M_I*om*kap;
    kcpl1 = num2*Kabbala(string("\\sqr(M_{")+flh.TexName()+string("})"),sqr(flh.Yuk()));
    //kcpl0=Kabbala("",300);
    //kcpl1=Kabbala("",301);
    
    v[vanz].cpl[0]  = kcpl0.Value();
    v[vanz].cpl[1]  = kcpl1.Value();
    v[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
    v[vanz].cpl[2]  = 0.;v[vanz].cpl[3]  = 0.;

    v[vanz].ncf   = 1;
    v[vanz].Color = new Color_Function; 
    
    v[vanz].Color->type   = cf::None;     
    
    v[vanz].nlf     = 1;
    v[vanz].Lorentz = new Lorentz_Function; 
    
    v[vanz].Lorentz->type = lf::SSGS;     
    v[vanz].Lorentz->SetParticleArg(0,2);     	

    v[vanz].on      = 1;
    vanz++;
  }
}


void Model_EW_Grav::c_SSST(Single_Vertex* v,int& vanz)
{
  Kabbala kcpl0,kcpl1,yuk;
  Kabbala num3  = Kabbala(string("3"),3.);

  Flavour flgraviton(kf::graviton);
  Flavour flgs(kf::gscalar);
  Flavour flh = Flavour(kf::h);

 if (flh.IsOn()&&flgraviton.IsOn()) {  
    v[vanz].nleg     = 4;
    v[vanz].in[0] = flh;
    v[vanz].in[1] = flh;
    v[vanz].in[2] = flh;
    v[vanz].in[3] = flgraviton;

    yuk   = Kabbala(string("M_{")+flh.TexName()+string("}"),flh.Yuk());
    kcpl0 = -M_I*kap*yuk*yuk*(num3/vev);
    kcpl1 = kcpl0;
    //kcpl0=Kabbala("",320);
    //kcpl1=Kabbala("",321);
    
    v[vanz].cpl[0]  = kcpl0.Value();
    v[vanz].cpl[1]  = kcpl1.Value();
    v[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
    v[vanz].cpl[2]  = 0.;v[vanz].cpl[3]  = 0.;

    v[vanz].ncf   = 1;
    v[vanz].Color = new Color_Function; 
    
    v[vanz].Color->type   = cf::None;     
    
    v[vanz].nlf     = 1;
    v[vanz].Lorentz = new Lorentz_Function; 
    
    v[vanz].Lorentz->type = lf::SSST;     
    v[vanz].Lorentz->SetParticleArg(3);     	

    v[vanz].on      = 1;
    vanz++;
  }

 if (flh.IsOn()&&flgs.IsOn()) {  
    v[vanz].nleg     = 4;
    v[vanz].in[0] = flh;
    v[vanz].in[1] = flh;
    v[vanz].in[2] = flh;
    v[vanz].in[3] = flgs;

    yuk   = Kabbala(string("M_{")+flh.TexName()+string("}"),flh.Yuk());
    kcpl0 = -M_I*num2*om*kap*yuk*yuk*(num3/vev);
    kcpl1 = kcpl0;
    //kcpl0=Kabbala("",300);
    //kcpl1=Kabbala("",0);
    
    v[vanz].cpl[0]  = kcpl0.Value();
    v[vanz].cpl[1]  = kcpl1.Value();
    v[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
    v[vanz].cpl[2]  = 0.;v[vanz].cpl[3]  = 0.;

    v[vanz].ncf   = 1;
    v[vanz].Color = new Color_Function; 
    
    v[vanz].Color->type   = cf::None;     
    
    v[vanz].nlf     = 1;
    v[vanz].Lorentz = new Lorentz_Function; 
    
    v[vanz].Lorentz->type = lf::SSSS;     

    v[vanz].on      = 1;
    vanz++;
  }
}


inline double Model_EW_Grav::Aqed(double t) {return aqed->Aqed(t);}
inline double Model_EW_Grav::Aqed()         {return aqed->AqedFixed();}

/*
inline double Model_EW::Aqed(double t) {return (*APHYTOOLS::aqed)(t);}

inline double Model_EW::Aqed() {return APHYTOOLS::aqed->Aqed(sqr(rpa.gen.Ecms()));}
*/

inline Complex Model_EW_Grav::CKM(short int i,short int j)   {return CplEW.CKM(i,j);}
inline Kabbala Model_EW_Grav::K_CKM(short int i,short int j)       
{   
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("V_{")+string(hi)+string(hj)+string("}"),CplEW.CKM(i,j));
} 
  
inline Kabbala Model_EW_Grav::conj_K_CKM(short int i,short int j)       
{   
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("V^\\ti_{")+string(hi)+string(hj)+string("}"),conj(CplEW.CKM(i,j)));
} 
 



  


