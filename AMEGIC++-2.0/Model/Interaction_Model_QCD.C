#include "Interaction_Model_QCD.H"
#include "MathTools.H"
#include "Message.H"
#include "Run_Parameter.H"

using namespace AMEGIC;
using namespace ATOOLS;
using namespace std;

Interaction_Model_QCD::Interaction_Model_QCD(MODEL::Model_Base * _model,
					     std::string _cplscheme,std::string _yukscheme) :
  Interaction_Model_Base(_model,_cplscheme,_yukscheme)
{ 
  double Ecms2 = sqr(rpa.gen.Ecms());

  g3  = Kabbala(string("g_3"),sqrt(4.*M_PI*ScalarFunction(std::string("alpha_S"),Ecms2)));
  PL  = Kabbala(string("P_L"),1.);
  PR  = Kabbala(string("P_R"),1.);
  M_I = Kabbala(string("i"),Complex(0.,1.)); 
}

void Interaction_Model_QCD::c_FFV(Single_Vertex * vertex,int & vanz)
{
  Kabbala kcpl0 = -g3*M_I;
  Kabbala kcpl1 = kcpl0;

  for (short int i=1;i<=6;i++) {
    Flavour flav = Flavour(kf::code(i));
    if (flav.Strong() && flav.IsOn() && Flavour(kf::gluon).IsOn()) { 
      vertex[vanz].in[0]         = flav;
      vertex[vanz].in[1]         = Flavour(kf::gluon);
      vertex[vanz].in[2]         = flav;

      vertex[vanz].cpl[0]        = kcpl0.Value();
      vertex[vanz].cpl[1]        = kcpl1.Value();
      vertex[vanz].cpl[2]        = 0.;
      vertex[vanz].cpl[3]        = 0.;
      vertex[vanz].Str           = (kcpl0*PR+kcpl1*PL).String();      

      vertex[vanz].ncf           = 1;
      vertex[vanz].Color         = new Color_Function(cf::T,1,2,0,'1','2','0');     

      vertex[vanz].nlf           = 1;
      vertex[vanz].Lorentz       = new Lorentz_Function(lf::Gamma);     
      vertex[vanz].Lorentz->SetParticleArg(1);     
                  
      vertex[vanz].on            = 1;
      vanz++;
    } 
  }
}

void Interaction_Model_QCD::c_VVV(Single_Vertex * vertex,int& vanz)
{
  Kabbala kcpl0 = -g3;
  Kabbala kcpl1 = kcpl0; 
  
  if (Flavour(kf::gluon).IsOn()) {

  for (short int i=0;i<3;i++) vertex[vanz].in[i] = Flavour(kf::gluon);

  vertex[vanz].cpl[0]        = kcpl0.Value();
  vertex[vanz].cpl[1]        = kcpl1.Value();
  vertex[vanz].cpl[2]        = 0.;
  vertex[vanz].cpl[3]        = 0.;
  vertex[vanz].Str           = (kcpl0*PR+kcpl1*PL).String();

  vertex[vanz].ncf           = 1;
  vertex[vanz].Color         = new Color_Function(cf::F);     
  vertex[vanz].Color->SetParticleArg(0,2,1);     
  vertex[vanz].Color->SetStringArg('0','2','1');     

  vertex[vanz].nlf           = 1;
  vertex[vanz].Lorentz       = new Lorentz_Function(lf::Gauge3);     
  vertex[vanz].Lorentz->SetParticleArg(0,1,2);     

  vertex[vanz].on            = 1;
  vanz++;

  }
}

void Interaction_Model_QCD::c_VVVV(Single_Vertex * vertex,int& vanz)
{
  Kabbala kcpl0 = -M_I*g3*g3; 
  Kabbala kcpl1 = kcpl0; 
  
  if (Flavour(kf::gluon).IsOn()) { 

  for (short int i=0;i<4;i++) vertex[vanz].in[i] = Flavour(kf::gluon);

  vertex[vanz].nleg            = 4;
  vertex[vanz].cpl[0]          = kcpl0.Value();
  vertex[vanz].cpl[1]          = kcpl1.Value();
  vertex[vanz].cpl[2]          = 0.;
  vertex[vanz].cpl[3]          = 0.;
  vertex[vanz].Str             = (kcpl0*PR+kcpl1*PL).String();
  
  vertex[vanz].ncf             = 3;
  vertex[vanz].nlf             = 3;
  
  vertex[vanz].Color           = new Color_Function[3];
  vertex[vanz].Lorentz         = new Lorentz_Function[3]; 

  vertex[vanz].Color[0]        = Color_Function(cf::F,0,2,4,'0','2','4',
				   new Color_Function(cf::F,1,3,4,'1','3','4'));
  vertex[vanz].Lorentz[0]      = Lorentz_Function(lf::Gluon4);
  vertex[vanz].Lorentz[0].SetParticleArg(0,1,2,3);     

  vertex[vanz].Color[1]        = Color_Function(cf::F,0,3,4,'0','3','4',
				   new Color_Function(cf::F,1,2,4,'1','2','4'));
  vertex[vanz].Lorentz[1]      = Lorentz_Function(lf::Gluon4);
  vertex[vanz].Lorentz[1].SetParticleArg(0,1,3,2);     

  vertex[vanz].Color[2]        = Color_Function(cf::F,0,1,4,'0','1','4',
				   new Color_Function(cf::F,3,2,4,'3','2','4')); 
  vertex[vanz].Lorentz[2]      = Lorentz_Function(lf::Gluon4);     
  vertex[vanz].Lorentz[2].SetParticleArg(0,3,1,2);     
  
  vertex[vanz].on              = 1;
  vanz++;
  }  
}

