#include "Model_QCD_Grav.H"
#include "MathTools.H"
#include "Run_Parameter.H"
#include "Running_AlphaS.H"

using namespace AMEGIC;
using namespace AORGTOOLS;
using namespace AMATOOLS;
using namespace APHYTOOLS;
using namespace std;

void Model_QCD_Grav::Init()
{
  // Spectrum
  //SpQCD.FillYukawas();
  //Couplings
  //Cpl.Init();
  CplLED.Init();
  g3 = Kabbala(string("g_3"),sqrt(4.*M_PI*Aqcd()));
  PL = Kabbala(string("P_L"),1.);
  PR = Kabbala(string("P_R"),1.);
  M_I = Kabbala(string("i"),Complex(0.,1.)); 
  kap   = Kabbala(string("kappa"),CplLED.Kappa());
  om    = Kabbala(string("omega"),CplLED.Omega());
  num2  = Kabbala(string("2"),2.);
  num15 = Kabbala(string("1.5"),1.5);
}

void Model_QCD_Grav::c_VVT(Single_Vertex* v,int& vanz)
{
  Kabbala kcpl0,kcpl1; 
  
  Flavour flgraviton(kf::graviton);
  Flavour flgs(kf::gscalar);
  Flavour flav(kf::gluon);
  if(flgraviton.IsOn()){

    v[vanz].in[0] = flav;
    v[vanz].in[1] = flgraviton;
    v[vanz].in[2] = flav;
    
    kcpl0 = -M_I*kap;
    kcpl1 = Kabbala();
    
    v[vanz].cpl[0]  = kcpl0.Value();
    v[vanz].cpl[1]  = kcpl1.Value();
    v[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
    v[vanz].cpl[2]  = 0.;v[vanz].cpl[3]  = 0.;

    v[vanz].ncf   = 1;
    v[vanz].Color = new Color_Function; 

    v[vanz].Color->type       = cf::G;//GD;     
    v[vanz].Color->SetParticleArg(0,2);     
    v[vanz].Color->SetStringArg('0','2');     
    
    v[vanz].nlf     = 1;
    v[vanz].Lorentz = new Lorentz_Function; 
  
    v[vanz].Lorentz->type = lf::VVT;     
    v[vanz].Lorentz->SetParticleArg(0,2,1);     

    v[vanz].on      = 1;
    vanz++;
  }

  // gluon gscalar gluon
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
    
    v[vanz].Color->type       = cf::G;//GD;     
    v[vanz].Color->SetParticleArg(0,2);     
    v[vanz].Color->SetStringArg('0','2');     
    
    v[vanz].nlf     = 1;
    v[vanz].Lorentz = new Lorentz_Function; 

    v[vanz].Lorentz->type = lf::VVGS;     
    v[vanz].Lorentz->SetParticleArg(0,2,1);     

    v[vanz].on      = 1;
    vanz++;
  }
}

void Model_QCD_Grav::c_FFVT(Single_Vertex* v,int& vanz)
{
  
  Flavour flgraviton(kf::graviton);
  Flavour flgs(kf::gscalar);
  Kabbala kcpl0,kcpl1;
  
  for (short int i=1;i<7;i++) {
    Flavour flav = Flavour(kf::code(i));
    if (flav.Strong() && flav.IsOn()) {
      if(flgraviton.IsOn()) { 
	v[vanz].nleg    = 4;
	v[vanz].in[0] = flav;
	v[vanz].in[1] = Flavour(kf::gluon);
	v[vanz].in[2] = flav;
	v[vanz].in[3] = flgraviton;
	
	kcpl0 = g3*M_I*kap/num2;
	kcpl1 = kcpl0;
	
	v[vanz].cpl[0]  = kcpl0.Value();
	v[vanz].cpl[1]  = kcpl1.Value();
	v[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	
	v[vanz].cpl[2]  = 0.;
	v[vanz].cpl[3]  = 0.;
	
	v[vanz].ncf   = 1;
	v[vanz].Color = new Color_Function; 
	
	v[vanz].Color->type       = cf::T;     
	//v[vanz].Color->SetParticleArg(1,0,2);     
	//v[vanz].Color->SetStringArg('1','0','2');     
	// To be tested !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	v[vanz].Color->SetParticleArg(1,2,0);     
	v[vanz].Color->SetStringArg('1','2','0');     
	
	//AORGTOOLS::msg.Out()<<"Test: "<<v[vanz].Color->String()<<endl;
	
	v[vanz].nlf     = 1;
	v[vanz].Lorentz = new Lorentz_Function; 

	v[vanz].Lorentz->type       = lf::FFVT;     
	v[vanz].Lorentz->SetParticleArg(1,3);     
	
	v[vanz].on      = 1;
	vanz++;
      }
      if(flgs.IsOn()){
	v[vanz].nleg    = 4;
	v[vanz].in[0] = flav;
	v[vanz].in[1] = Flavour(kf::gluon);
	v[vanz].in[2] = flav;
	v[vanz].in[3] = flgs;
	
	kcpl0 = -g3*M_I*om*kap*num15;
	kcpl1 = kcpl0;
	
	v[vanz].cpl[0]  = kcpl0.Value();
	v[vanz].cpl[1]  = kcpl1.Value();
	v[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	
	v[vanz].cpl[2]  = 0.;
	v[vanz].cpl[3]  = 0.;
	
	v[vanz].ncf   = 1;
	v[vanz].Color = new Color_Function; 
	
	v[vanz].Color->type       = cf::T;     
	//v[vanz].Color->SetParticleArg(1,0,2);     
	//v[vanz].Color->SetStringArg('1','0','2');     
	// To be tested !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	v[vanz].Color->SetParticleArg(1,2,0);     
	v[vanz].Color->SetStringArg('1','2','0');     
	
	//AORGTOOLS::msg.Out()<<"Test: "<<v[vanz].Color->String()<<endl;
	
	v[vanz].nlf     = 1;
	v[vanz].Lorentz = new Lorentz_Function; 

	v[vanz].Lorentz->type       = lf::FFVGS;     
	v[vanz].Lorentz->SetParticleArg(1,3);     
	
	v[vanz].on      = 1;
	vanz++;
      }
    } 
  }
}

void Model_QCD_Grav::c_VVVT(Single_Vertex* v,int& vanz)
{
  Kabbala kcpl0,kcpl1; 
  
  Flavour flgraviton(kf::graviton);
  if(!flgraviton.IsOn())return;

  v[vanz].nleg    = 4;
  for (short int i=0;i<3;i++)
    v[vanz].in[i] = Flavour(kf::gluon);
  v[vanz].in[3] = flgraviton;
  

  kcpl0 = g3*kap; 
  kcpl1 = kcpl0; 

  v[vanz].cpl[0]  = kcpl0.Value();
  v[vanz].cpl[1]  = kcpl1.Value();
  v[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
  v[vanz].cpl[2]  = 0.;v[vanz].cpl[3]  = 0.;

  v[vanz].ncf   = 1;
  v[vanz].Color = new Color_Function; 

  v[vanz].Color->type       = cf::F;     
  v[vanz].Color->SetParticleArg(0,2,1);     
  v[vanz].Color->SetStringArg('0','2','1');     

  v[vanz].nlf     = 1;
  v[vanz].Lorentz = new Lorentz_Function; 

  v[vanz].Lorentz->type = lf::VVVT;     
  v[vanz].Lorentz->SetParticleArg(0,1,2,3);     

  v[vanz].on      = 1;
  vanz++;
}


inline double Model_QCD_Grav::Aqcd(double t)    {return APHYTOOLS::as->AlphaS(t);} // switch dependent running
inline double Model_QCD_Grav::Aqcd()            {return APHYTOOLS::as->AsFixed();} // alpha_S _eff (read in)

/*
inline double Model_QCD::Aqcd(double t) {return (*as)(t);}
inline double Model_QCD::Aqcd()         {return (*as)(sqr(rpa.gen.Ecms()));}
*/


















