#include "Model_QCD.H"
#include "MathTools.H"
#include "Run_Parameter.H"
#include "Message.H"
#include "Running_AlphaS.H"

using namespace AMEGIC;
using namespace AORGTOOLS;
using namespace APHYTOOLS;
using namespace AMATOOLS;
using namespace std;

void Model_QCD::Init()
{
  // Spectrum
  //  SpQCD.Fill_Masses();
  //Couplings
  //Cpl.Init();
  g3 = Kabbala(string("g_3"),sqrt(4.*M_PI*Aqcd()));
  PL = Kabbala(string("P_L"),1.);
  PR = Kabbala(string("P_R"),1.);
  M_I = Kabbala(string("i"),Complex(0.,1.)); 
}

void Model_QCD::c_FFV(Single_Vertex* v,int& vanz)
{
  
  Kabbala kcpl0,kcpl1;
  
  for (short int i=1;i<=6;i++) {
    Flavour flav = Flavour(kf::code(i));
    if (flav.Strong() && flav.IsOn() && Flavour(kf::gluon).IsOn()) { 
      v[vanz].in[0] = flav;
      v[vanz].in[1] = Flavour(kf::gluon);
      v[vanz].in[2] = flav;
      
      kcpl0 = -g3*M_I;
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
      

      v[vanz].nlf     = 1;
      v[vanz].Lorentz = new Lorentz_Function; 

      v[vanz].Lorentz->type       = lf::Gamma;     
      v[vanz].Lorentz->SetParticleArg(1);     
                  
      v[vanz].on      = 1;
      vanz++;
      
    } 
  }
}

void Model_QCD::c_VVV(Single_Vertex* v,int& vanz)
{
  Kabbala kcpl0,kcpl1; 
 
  if (Flavour(kf::gluon).IsOn()) { 
  
  for (short int i=0;i<3;i++)
    v[vanz].in[i] = Flavour(kf::gluon);
  
  kcpl0 = -g3; 
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

  v[vanz].Lorentz->type = lf::Gauge3;     
  v[vanz].Lorentz->SetParticleArg(0,1,2);     

  v[vanz].on      = 1;
  vanz++;
  }
}

void Model_QCD::c_VVVV(Single_Vertex* v,int& vanz)
{
  Kabbala kcpl0,kcpl1; 

  if (Flavour(kf::gluon).IsOn()) { 
  
  for (short int i=0;i<4;i++)
    v[vanz].in[i] = Flavour(kf::gluon);
  
  kcpl0 = -M_I*g3*g3; 
  kcpl1 = kcpl0; 
  
  v[vanz].nleg    = 4;
  v[vanz].cpl[0]  = kcpl0.Value();
  v[vanz].cpl[1]  = kcpl1.Value();
  v[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
  v[vanz].cpl[2]  = 0.;v[vanz].cpl[3]  = 0.;
  
  v[vanz].ncf   = 3;
  v[vanz].nlf   = 3;
  
  v[vanz].Color   = new Color_Function[3];
  v[vanz].Lorentz = new Lorentz_Function[3]; 

  v[vanz].Color[0]        = Color_Function(cf::F,0,2,4,'0','2','4');
  (v[vanz].Color[0]).Next = new Color_Function(cf::F,1,3,4,'1','3','4');
  v[vanz].Lorentz[0].type = lf::Gluon4;      
  v[vanz].Lorentz[0].SetParticleArg(0,1,2,3);     
    
  v[vanz].Color[1]        = Color_Function(cf::F,0,3,4,'0','3','4');
  (v[vanz].Color[1]).Next = new Color_Function(cf::F,1,2,4,'1','2','4');
  v[vanz].Lorentz[1].type = lf::Gluon4;      
  v[vanz].Lorentz[1].SetParticleArg(0,1,3,2);     

  v[vanz].Color[2]        = Color_Function(cf::F,0,1,4,'0','1','4');
  (v[vanz].Color[2]).Next = new Color_Function(cf::F,3,2,4,'3','2','4'); 
  v[vanz].Lorentz[2].type = lf::Gluon4;     
  v[vanz].Lorentz[2].SetParticleArg(0,3,1,2);     
  
  v[vanz].on      = 1;
  vanz++;
  }  
}


inline double Model_QCD::Aqcd(double t) { return as->AlphaS(t); } // switch dependent running
inline double Model_QCD::Aqcd()         { return as->AsFixed(); } // alpha_S _eff (read in)

