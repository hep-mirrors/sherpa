#include "Model_EE_QCD.H"
#include "Run_Parameter.H"
#include "Running_AlphaS.H"
#include "Running_AlphaQED.H"
#include "MathTools.H"

using namespace AMEGIC;
using namespace AORGTOOLS;
using namespace APHYTOOLS;
using namespace AMATOOLS;
using namespace std;

void Model_EE_QCD::Init()
{
  //QCD part
  moqcd.Init();
  //Couplings
  Cpl.Init();
  
  g1 = Kabbala(string("g_1"),sqrt(4.*M_PI*Aqed()));
  g2 = Kabbala(string("g_1/sin\\theta_W"), g1.Value()/SinTW());
  PL = Kabbala(string("P_L"),1.);
  PR = Kabbala(string("P_R"),1.);
  M_I = Kabbala(string("i"),Complex(0.,1.));
}

void Model_EE_QCD::c_FFV(Single_Vertex* v,int& vanz)
  {
  moqcd.c_FFV(v,vanz);
  Flavour flphoton(kf::photon);
  Flavour flZ(kf::Z);
    
  Kabbala kcpl0,kcpl1;
  
  // only electrons
  for (short int i=1;i<12;i++) {
    if (i==7) i=11;
    
    Flavour flav = Flavour(kf::code(i));
    Kabbala charge = Kabbala(string("Q_{")+ string(flav.texname())+string("}"),flav.charge());
    Kabbala isoweak = Kabbala(string("T_{")+ string(flav.texname())+string("}"),flav.isoweak());
    

    if (flav.ison()) { 
      //photon
      if (flphoton.ison()) {
	v[vanz].in[0] = flav;
	v[vanz].in[2] = flav;
	v[vanz].in[1] = Flavour(kf::photon);
	
	kcpl0 = -g1*M_I*charge;
	kcpl1 = kcpl0;
	
	v[vanz].cpl[0]  = kcpl0.Value();
	v[vanz].cpl[1]  = kcpl1.Value();
	v[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	
	v[vanz].cpl[2]  = 0.;
	v[vanz].cpl[3]  = 0.;

	v[vanz].ncf   = 1;
	v[vanz].Color = new Color_Function; 
	
	if (flav.strong()) {
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
      //Z
      if (flZ.ison()) {
	v[vanz].in[0] = flav;
	v[vanz].in[2] = flav;
	v[vanz].in[1] = Flavour(kf::Z);
	
	kcpl0 = -M_I*(isoweak-charge*K_sinTW()*K_sinTW())*g2/K_cosTW();
	kcpl1 = M_I*charge*K_sinTW()*K_sinTW()*g2/K_cosTW();
	
	v[vanz].cpl[0]  = kcpl0.Value();
	v[vanz].cpl[1]  = kcpl1.Value();
	v[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();

	v[vanz].cpl[2]  = 0.;
	v[vanz].cpl[3]  = 0.;
	
	v[vanz].ncf   = 1;
	v[vanz].Color = new Color_Function; 
	
	if (flav.strong()) {
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

void Model_EE_QCD::c_VVV(Single_Vertex* v,int& vanz)
{
  moqcd.c_VVV(v,vanz);
}
void Model_EE_QCD::c_VVVV(Single_Vertex* v,int& vanz)
{
  moqcd.c_VVVV(v,vanz);
}

inline double Model_EE_QCD::Aqcd(double t) {return as->AlphaS(t);} // switch dependent running
inline double Model_EE_QCD::Aqcd()         {return as->AsFixed();} // alpha_S _eff (read in)
inline double Model_EE_QCD::Aqed(double t) {return aqed->Aqed(t);}
inline double Model_EE_QCD::Aqed()         {return aqed->AqedFixed();}

/*
inline double Model_EE_QCD::Aqcd(double t) {return (*as)(t);}
inline double Model_EE_QCD::Aqcd()         {return (*as)(sqr(rpa.gen.Ecms()));}
inline double Model_EE_QCD::Aqed(double t) {return (*aqed)(t);}
inline double Model_EE_QCD::Aqed()         {return (*aqed)(sqr(rpa.gen.Ecms()));}
*/














