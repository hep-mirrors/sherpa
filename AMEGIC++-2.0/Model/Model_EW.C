#include <stdio.h>

#include "Model_EW.H"
#include "Running_AlphaQED.H"
#include "Run_Parameter.H"
#include "MathTools.H"
#include "Vector.H"

using namespace AMEGIC;
using namespace AORGTOOLS;
using namespace APHYTOOLS;
using namespace AMATOOLS;
using namespace std;

void Model_EW::Init()
{
  // Initialize Yukawa masses
  SpEW.FillYukawas();
  //Couplings
  CplEW.Init();

  g1    = Kabbala(string("g_1"),sqrt(4.*M_PI*Aqed()));
  g2    = Kabbala(string("g_1/sin\\theta_W"), g1.Value()/SinTW());
  PL    = Kabbala(string("P_L"),1.);
  PR    = Kabbala(string("P_R"),1.);
  M_I   = Kabbala(string("i"),Complex(0.,1.));
  root2 = Kabbala(string("\\sqrt{2}"),sqrt(2.));
  //vev   = Kabbala(string("v_{EW}"),CplEW.VEV());
  vev   = Kabbala(string("v_{EW}"),2.*SinTW()*CosTW()*Flavour(kf::Z).Mass()/g1.Value());
}

void Model_EW::c_FFV(Single_Vertex* v,int& vanz)
{
  Flavour flphoton(kf::photon);
  Flavour flZ(kf::Z);
  Flavour flW(kf::W);
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
	    //photon
	    if (flphoton.IsOn()) {
	      
	      Kabbala kcpl0,kcpl1;
	      kcpl0 = -g1*M_I*charge1;
	      kcpl1 = kcpl0;
	      
	      if (!AMATOOLS::IsZero(kcpl0.Value())) {
		v[vanz].in[0] = flav1;
		v[vanz].in[1] = Flavour(kf::photon);
		v[vanz].in[2] = flav2;
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
		
		v[vanz].Lorentz->type       = lf::Gamma;     
		v[vanz].Lorentz->SetParticleArg(1);     

		v[vanz].on      = 1;
		vanz++;
	      }
	    }
	    //Z
	    if (flZ.IsOn()) {
	      
	      Kabbala kcpl0,kcpl1;
	      kcpl0 = M_I/K_cosTW()*charge1*K_sinTW()*K_sinTW()*g2;
	      kcpl1 = -M_I/K_cosTW()*(isoweak1-charge1*K_sinTW()*K_sinTW())*g2;
	     
	      v[vanz].in[0] = flav1;
	      v[vanz].in[1] = Flavour(kf::Z);
	      v[vanz].in[2] = flav2;
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
	      
	      v[vanz].Lorentz->type       = lf::Gamma;     
	      v[vanz].Lorentz->SetParticleArg(1);     

	      v[vanz].on     = 1;
	      vanz++;
	    }
	  }
	  //W
	  if (flW.IsOn()) {
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
		kcpl1 = -M_I/root2*g2;
	      if (i<7 && j<7) {
		if (flav1.IsDowntype())
		kcpl1 = -M_I/root2*g2*K_CKM((i-1)/2,j/2-1);
		else 	    
		kcpl1 = -M_I/root2*g2*K_CKM(i/2-1,(j-1)/2);		
	      }
	      if (!AMATOOLS::IsZero(kcpl1.Value())) {
		v[vanz].in[1] = Flavour(kf::W);
		if (flav1.IsDowntype()) {
		  v[vanz].in[0] = flav1;
		  v[vanz].in[2] = flav2;
		}
		else {
		  v[vanz].in[0] = flav2;
		  v[vanz].in[2] = flav1;
		}
		
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
		
		v[vanz].Lorentz->type       = lf::Gamma;     
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

void Model_EW::c_VVV(Single_Vertex* v,int& vanz)
{
  Flavour flav(kf::W);
  Kabbala kcpl0,kcpl1,kcpl0_1,kcpl1_1,charge;
  charge = Kabbala(string("Q_{")+flav.TexName()+string("}"),flav.Charge());

  if (!flav.IsOn()) return;

  // photon WW
  v[vanz].in[0] = flav;
  v[vanz].in[1] = Flavour(kf::photon);
  v[vanz].in[2] = flav;

  kcpl0 = M_I*g1*charge;
  kcpl1 = kcpl0;

  v[vanz].cpl[0]  = kcpl0.Value();
  v[vanz].cpl[1]  = v[vanz].cpl[0];
  v[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
  v[vanz].cpl[2]  = 0.;v[vanz].cpl[3]  = 0.;

  v[vanz].ncf   = 1;
  v[vanz].Color = new Color_Function; 

  v[vanz].Color->type       = cf::None;     

  v[vanz].nlf     = 1;
  v[vanz].Lorentz = new Lorentz_Function; 

  v[vanz].Lorentz->type = lf::Gauge3;     
  v[vanz].Lorentz->SetParticleArg(0,1,2);     

  v[vanz].on      = 1;
  vanz++;

  // ZWW
  v[vanz].in[0] = flav;
  v[vanz].in[1] = Flavour(kf::Z);
  v[vanz].in[2] = flav;

  kcpl0 = M_I*g2*charge*K_cosTW();
  kcpl1 = kcpl0;
  
  v[vanz].cpl[0]  = kcpl0.Value();
  v[vanz].cpl[1]  = kcpl1.Value();
  v[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
  v[vanz].cpl[2]  = 0.;v[vanz].cpl[3]  = 0.;

  v[vanz].ncf   = 1;
  v[vanz].Color = new Color_Function; 

  v[vanz].Color->type       = cf::None;     

  v[vanz].nlf     = 1;
  v[vanz].Lorentz = new Lorentz_Function; 

  v[vanz].Lorentz->type = lf::Gauge3;     
  v[vanz].Lorentz->SetParticleArg(0,1,2);     

  v[vanz].on      = 1;
  vanz++;
}


void Model_EW::c_FFS(Single_Vertex* v,int& vanz)
{
  Flavour flh(kf::h);
  Kabbala kcpl0,kcpl1,M_h;
  if (!flh.IsOn()) return;

  for (short int i=1;i<17;i++) {
    if (i==7) i=11;
    Flavour flav = Flavour(kf::code(i));
    if (flav.IsOn() && flav.IsFermion() && (flav.Yuk() > 0.05)) {
      
      M_h = Kabbala(string("M_{")+flav.TexName()+string("}"),flav.Yuk());

      kcpl0 = -M_I*M_h/vev;
      kcpl1 = kcpl0;
      
      if (!AMATOOLS::IsZero(kcpl0.Value())) {
	v[vanz].in[0] = flav;
	v[vanz].in[1] = flh;
	v[vanz].in[2] = flav;

	v[vanz].cpl[0]  = kcpl0.Value();
	v[vanz].cpl[1]  = kcpl1.Value();
	v[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	v[vanz].cpl[2]  = 0.;v[vanz].cpl[3]  = 0.;

	v[vanz].ncf   = 1;
	v[vanz].Color = new Color_Function; 
	
	if (flav.Strong()) {
	  v[vanz].Color->type       = cf::D;     
	  v[vanz].Color->SetParticleArg(0,2);     
	  v[vanz].Color->SetStringArg('0','2');     
	}
	else v[vanz].Color->type = cf::None; 
	
	v[vanz].nlf     = 1;
	v[vanz].Lorentz = new Lorentz_Function; 
	
	v[vanz].Lorentz->type = lf::FFS;     

	v[vanz].on      = 1;
	vanz++;
      }
    }
  }
}


void Model_EW::c_VVS(Single_Vertex* v,int& vanz)
{
  Flavour flh(kf::h);
  Kabbala kcpl0,kcpl1;  
  Kabbala num_2 = Kabbala(string("2"),2.);  
 
 if (!flh.IsOn()) return;
  
  Flavour flav(kf::W);
  // W h W
  if (flav.IsOn()) {
    v[vanz].in[0] = flav;
    v[vanz].in[1] = flh;
    v[vanz].in[2] = flav;
    
    kcpl0 = M_I*g2*flav.Yuk()/vev;
    kcpl1 = kcpl0;
    
    v[vanz].cpl[0]  = kcpl0.Value();
    v[vanz].cpl[1]  = v[vanz].cpl[0];
    v[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
    v[vanz].cpl[2]  = 0.;v[vanz].cpl[3]  = 0.;

    v[vanz].ncf   = 1;
    v[vanz].Color = new Color_Function; 
    
    v[vanz].Color->type       = cf::None;     

    v[vanz].nlf     = 1;
    v[vanz].Lorentz = new Lorentz_Function; 

    v[vanz].Lorentz->type = lf::Gab;     
    v[vanz].Lorentz->SetParticleArg(0,2);     

    v[vanz].on      = 1;
    vanz++;
  }

  flav = Flavour(kf::Z);
  // Z h Z
  if (flav.IsOn()) {
    v[vanz].in[0] = flav;
    v[vanz].in[1] = flh;
    v[vanz].in[2] = flav;
    
    kcpl0 = M_I*g2*flav.Yuk()/vev;
    kcpl1 = kcpl0;

    v[vanz].cpl[0]  = kcpl0.Value();
    v[vanz].cpl[1]  = kcpl1.Value();
    v[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
    v[vanz].cpl[2]  = 0.;v[vanz].cpl[3]  = 0.;

    v[vanz].ncf   = 1;
    v[vanz].Color = new Color_Function; 
    
    v[vanz].Color->type       = cf::None;     

    v[vanz].nlf     = 1;
    v[vanz].Lorentz = new Lorentz_Function; 

    v[vanz].Lorentz->type = lf::Gab;     
    v[vanz].Lorentz->SetParticleArg(0,2);     

    v[vanz].on      = 1;
    vanz++;
  }
}


void Model_EW::c_SSS(Single_Vertex* v,int& vanz)
{
  Flavour flh = Flavour(kf::h);
  Kabbala kcpl0,kcpl1,yuk;  
  Kabbala num_3 = Kabbala(string("3"),3.);  

 if (flh.IsOn()) {  
    v[vanz].in[0] = flh;
    v[vanz].in[1] = flh;
    v[vanz].in[2] = flh;

    yuk   = Kabbala(string("M_{")+flh.TexName()+string("}"),flh.Yuk());
    kcpl0 = -M_I*yuk*yuk*(num_3/vev);
    kcpl1 = kcpl0;
    
    v[vanz].cpl[0]  = kcpl0.Value();
    v[vanz].cpl[1]  = kcpl1.Value();
    v[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
    v[vanz].cpl[2]  = 0.;v[vanz].cpl[3]  = 0.;

    v[vanz].ncf   = 1;
    v[vanz].Color = new Color_Function; 
    
    v[vanz].Color->type   = cf::None;     
    
    v[vanz].nlf     = 1;
    v[vanz].Lorentz = new Lorentz_Function; 
    
    v[vanz].Lorentz->type = lf::SSS;     

    v[vanz].on      = 1;
    vanz++;
  }
}

void Model_EW::c_SSSS(Single_Vertex* v,int& vanz)
{
  Flavour flh = Flavour(kf::h);
  Kabbala kcpl0,kcpl1,yuk;  
  Kabbala num_3 = Kabbala(string("3"),3.);  

  if (flh.IsOn()) {  
    v[vanz].in[0] = flh;
    v[vanz].in[1] = flh;
    v[vanz].in[2] = flh;
    v[vanz].in[3] = flh;

    v[vanz].nleg  = 4;  
    
    yuk   = Kabbala(string("M_{")+flh.TexName()+string("}"),flh.Yuk());
    kcpl0 = -M_I*yuk*yuk*(num_3/(vev*vev));
    kcpl1 = kcpl0;
    
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

void Model_EW::c_VVVV(Single_Vertex* v,int& vanz)
{

  Flavour flavW(kf::W);
  Flavour flavZ(kf::Z);
  Flavour flavP(kf::photon);
  Kabbala kcpl0,kcpl1;
  
  // Ph - W - W - Ph  
  if (flavW.IsOn() && flavP.IsOn()) {
  v[vanz].in[0] = flavP;
  v[vanz].in[1] = flavW;
  v[vanz].in[2] = flavW.Bar();
  v[vanz].in[3] = flavP;
  
  v[vanz].nleg     = 4;

  kcpl0 = -M_I*g1*g1;
  kcpl1 = kcpl0;
    
  v[vanz].cpl[0]  = kcpl0.Value();
  v[vanz].cpl[1]  = kcpl1.Value();
  v[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
  v[vanz].cpl[2]  = 0.;v[vanz].cpl[3]  = 0.;
  
  v[vanz].ncf   = 1;
  v[vanz].Color = new Color_Function; 
  
  v[vanz].Color->type   = cf::None;     
    
  v[vanz].nlf     = 1;
  v[vanz].Lorentz = new Lorentz_Function; 
  
  v[vanz].Lorentz->type = lf::Gauge4;     
  v[vanz].Lorentz->SetParticleArg(0,3,1,2);     

  v[vanz].on      = 1;
  vanz++;
  }

  // Ph - W - W - Z  
  if (flavW.IsOn() && flavP.IsOn() && flavZ.IsOn()) {
  v[vanz].in[0] = flavP;
  v[vanz].in[1] = flavW;
  v[vanz].in[2] = flavW.Bar();
  v[vanz].in[3] = flavZ;

  v[vanz].nleg     = 4;  

  kcpl0 = -M_I*g1*g1*K_cosTW()/K_sinTW();
  kcpl1 = kcpl0;
  
  v[vanz].cpl[0]  = kcpl0.Value();
  v[vanz].cpl[1]  = kcpl1.Value();
  v[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
  v[vanz].cpl[2]  = 0.;v[vanz].cpl[3]  = 0.;
  
  v[vanz].ncf   = 1;
  v[vanz].Color = new Color_Function; 
  
  v[vanz].Color->type   = cf::None;     
    
  v[vanz].nlf     = 1;
  v[vanz].Lorentz = new Lorentz_Function; 
  
  v[vanz].Lorentz->type = lf::Gauge4;     
  v[vanz].Lorentz->SetParticleArg(0,3,1,2);     

  v[vanz].on      = 1;
  //v[vanz].on      = 0;
  vanz++;
  }
  // Z - W - W - Z  
  if (flavW.IsOn() && flavZ.IsOn()) {
  v[vanz].in[0] = flavZ;
  v[vanz].in[1] = flavW;
  v[vanz].in[2] = flavW.Bar();
  v[vanz].in[3] = flavZ;
  
  v[vanz].nleg     = 4;

  kcpl0 = -M_I*g1*g1*K_cosTW()*K_cosTW()/(K_sinTW()*K_sinTW());
  kcpl1 = kcpl0;
  
  v[vanz].cpl[0]  = kcpl0.Value();
  v[vanz].cpl[1]  = kcpl1.Value();
  v[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
  v[vanz].cpl[2]  = 0.;v[vanz].cpl[3]  = 0.;
  
  v[vanz].ncf   = 1;
  v[vanz].Color = new Color_Function; 
  
  v[vanz].Color->type   = cf::None;     
    
  v[vanz].nlf     = 1;
  v[vanz].Lorentz = new Lorentz_Function; 
  
  v[vanz].Lorentz->type = lf::Gauge4;     
  v[vanz].Lorentz->SetParticleArg(0,3,1,2);     

  v[vanz].on      = 1;
  //v[vanz].on      = 0;
  vanz++;
  }
  
 // W - W - W - W  
  if (flavW.IsOn()) {
  v[vanz].in[0] = flavW.Bar();
  v[vanz].in[1] = flavW;
  v[vanz].in[2] = flavW.Bar();
  v[vanz].in[3] = flavW.Bar();
  
  v[vanz].nleg     = 4;

  kcpl0 = M_I*g1*g1/(K_sinTW()*K_sinTW());
  kcpl1 = kcpl0;
  
  v[vanz].cpl[0]  = kcpl0.Value();
  v[vanz].cpl[1]  = kcpl1.Value();
  v[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
  v[vanz].cpl[2]  = 0.;v[vanz].cpl[3]  = 0.;
  
  v[vanz].ncf   = 1;
  v[vanz].Color = new Color_Function; 
  
  v[vanz].Color->type   = cf::None;     
    
  v[vanz].nlf     = 1;
  v[vanz].Lorentz = new Lorentz_Function; 
  
  v[vanz].Lorentz->type = lf::Gauge4;     
  v[vanz].Lorentz->SetParticleArg(0,1,2,3);     

  v[vanz].on      = 1;
  vanz++;
  }
}

void Model_EW::c_SSVV(Single_Vertex* v,int& vanz)
{
  Kabbala num_2 = Kabbala(string("2"),2.);  

  Flavour flavW(kf::W);
  Flavour flavZ(kf::Z);
  Flavour flavh(kf::h);
  Kabbala kcpl0,kcpl1;
  
  // h - Z - Z - h  
  if (flavZ.IsOn() && flavh.IsOn()) {
    v[vanz].in[0] = flavZ;
    v[vanz].in[1] = flavh;
    v[vanz].in[2] = flavh;
    v[vanz].in[3] = flavZ;
    
    v[vanz].nleg     = 4;
    
    kcpl0 = M_I*g2*g2/(K_cosTW()*K_cosTW()*num_2);
    kcpl1 = kcpl0;
    
    v[vanz].cpl[0]  = kcpl0.Value();
    v[vanz].cpl[1]  = kcpl1.Value();
    v[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
    v[vanz].cpl[2]  = 0.;v[vanz].cpl[3]  = 0.;
    
    v[vanz].ncf   = 1;
    v[vanz].Color = new Color_Function; 
    
    v[vanz].Color->type   = cf::None;     
    
    v[vanz].nlf     = 1;
    v[vanz].Lorentz = new Lorentz_Function; 
    
    v[vanz].Lorentz->type = lf::VVSS;     
    v[vanz].Lorentz->SetParticleArg(0,3);     
    
    v[vanz].on      = 1;
    vanz++;
  }

  // h - W - W - h  
  if (flavW.IsOn() && flavh.IsOn()) {
    v[vanz].in[0] = flavW;
    v[vanz].in[1] = flavh;
    v[vanz].in[2] = flavh;
    v[vanz].in[3] = flavW;
    
    v[vanz].nleg     = 4;
    
    kcpl0 = M_I*g2*g2/num_2;
    kcpl1 = kcpl0;
    
    v[vanz].cpl[0]  = kcpl0.Value();
    v[vanz].cpl[1]  = kcpl1.Value();
    v[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
    v[vanz].cpl[2]  = 0.;v[vanz].cpl[3]  = 0.;
    
    v[vanz].ncf   = 1;
    v[vanz].Color = new Color_Function; 
    
    v[vanz].Color->type   = cf::None;     
    
    v[vanz].nlf     = 1;
    v[vanz].Lorentz = new Lorentz_Function; 
    
    v[vanz].Lorentz->type = lf::VVSS;     
    v[vanz].Lorentz->SetParticleArg(0,3);     
    
    v[vanz].on      = 1;
    vanz++;
  }


}

inline double Model_EW::Aqed(double t) {return aqed->Aqed(t);}
inline double Model_EW::Aqed()         {return aqed->AqedFixed();}

Complex Model_EW::CKM(short int i,short int j)   {return CplEW.CKM(i,j);}
Kabbala Model_EW::K_CKM(short int i,short int j)       
{   
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("V_{")+string(hi)+string(hj)+string("}"),CplEW.CKM(i,j));
} 
  
Kabbala Model_EW::conj_K_CKM(short int i,short int j)       
{   
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("V^\\ti_{")+string(hi)+string(hj)+string("}"),conj(CplEW.CKM(i,j)));
} 
 



  


