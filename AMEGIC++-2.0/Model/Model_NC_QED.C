#include "Model_NC_QED.H"
#include "Run_Parameter.H"
#include "MathTools.H"
#include "Running_AlphaQED.H"

using namespace AMEGIC;
using namespace AORGTOOLS;
using namespace APHYTOOLS;
using namespace AMATOOLS;
using namespace std;

void Model_NC_QED::Init()
{
  //Couplings
  Cpl.Init();

  //reading Theta-Matrix
  Data_Read dr(rpa.GetPath()+string("/")+rpa.me.ModelFile());

  double sinaNC = dr.GetValue<double>("sin(alpha_NC)");
  double cosaNC = sqrt(1.-sqr(sinaNC));
  double sinbNC = dr.GetValue<double>("sin(beta_NC)");
  double cosbNC = sqrt(1.-sqr(sinbNC));
  double singNC = dr.GetValue<double>("sin(gamma_NC)");
  double cosgNC = sqrt(1.-sqr(singNC));

  double lambda_NC = dr.GetValue<double>("Lambda_NC");
  
  Matrix<4>* Theta = rpa.me.GetTheta();

  for (short int i=0;i<4;i++) (*Theta)[i][i] = 0.;

  (*Theta)[0][1] = sinaNC*cosbNC;
  (*Theta)[0][2] = sinaNC*sinbNC;
  (*Theta)[0][3] = cosaNC;
  (*Theta)[1][2] = cosgNC;
  (*Theta)[1][3] = -singNC*sinbNC;
  (*Theta)[2][3] = -singNC*cosbNC;

  for (short int i=0;i<4;i++) {
    for (short int j=0;j<i;j++) (*Theta)[i][j] = -(*Theta)[j][i];
  }

  (*Theta) = 1./sqr(lambda_NC)*(*Theta);

  Theta->MatrixOut();

  g1   = Kabbala(string("g_1"),sqrt(4.*M_PI*Aqed()));
  M_I  = Kabbala(string("i"),Complex(0.,1.));

  PL = Kabbala(string("P_L"),1.);
  PR = Kabbala(string("P_R"),1.);

}

void Model_NC_QED::c_FFV(Single_Vertex* v,int& vanz)
{
  Flavour flphoton(kf::photon);
    
  Kabbala kcpl0,kcpl1;
  
  // only electrons
  for (short int i=11;i<16;i+=2) {
    
    Flavour flav    = Flavour(kf::code(i));
    Kabbala charge  = Kabbala(string("Q_{")+ string(flav.TexName())+string("}"),flav.Charge());
    Kabbala isoweak = Kabbala(string("T_{")+ string(flav.TexName())+string("}"),flav.IsoWeak());
    

    if (flav.IsOn()) { 
      //photon
      if (flphoton.IsOn()) {
	v[vanz].in[0] = flav;
	v[vanz].in[1] = Flavour(kf::photon);
	v[vanz].in[2] = flav;
	
	kcpl0 = -g1*M_I*charge;
	kcpl1 = kcpl0;
	
	v[vanz].cpl[0]  = kcpl0.Value();
	v[vanz].cpl[1]  = kcpl1.Value();
	v[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	
	v[vanz].cpl[2]  = 0.;
	v[vanz].cpl[3]  = 0.;

	v[vanz].ncf   = 1;
	v[vanz].Color = new Color_Function; 
	
	v[vanz].Color->type = cf::None; 
	
	v[vanz].nlf     = 1;
	v[vanz].Lorentz = new Lorentz_Function; 

	v[vanz].Lorentz->type = lf::Gamma_NC;     
	v[vanz].Lorentz->SetParticleArg(1);     

	v[vanz].on      = 1;
	vanz++;
      }
    }
  }
}


void Model_NC_QED::c_VVV(Single_Vertex* v,int& vanz)
{
  Kabbala kcpl0,kcpl1; 
  
  for (short int i=0;i<3;i++)
    v[vanz].in[i] = Flavour(kf::photon);
  
  kcpl0 = 2.*g1; 
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

  v[vanz].Lorentz->type = lf::Photon3_NC;     
  v[vanz].Lorentz->SetParticleArg(0,1,2);     

  v[vanz].on      = 1;
  vanz++;
}

void Model_NC_QED::c_VVVV(Single_Vertex* v,int& vanz)
{
  Kabbala kcpl0,kcpl1; 
  
  for (short int i=0;i<4;i++)
    v[vanz].in[i] = Flavour(kf::photon);
  
  kcpl0 = M_I*4*g1*g1; 
  kcpl1 = kcpl0; 
  
  v[vanz].nleg    = 4;
  v[vanz].cpl[0]  = kcpl0.Value();
  v[vanz].cpl[1]  = kcpl1.Value();
  v[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
  v[vanz].cpl[2]  = 0.;v[vanz].cpl[3]  = 0.;
  
  v[vanz].ncf   = 1;
  v[vanz].Color = new Color_Function; 
  v[vanz].Color->type       = cf::None;

  v[vanz].nlf     = 1;
  v[vanz].Lorentz = new Lorentz_Function; 
  v[vanz].Lorentz->type = lf::Photon4_NC;      
  v[vanz].Lorentz->SetParticleArg(0,1,2,3);     
    
  v[vanz].on      = 1;
  vanz++;  
}


inline double Model_NC_QED::Aqed(double t) {return aqed->Aqed(t);}
inline double Model_NC_QED::Aqed()         {return aqed->AqedFixed();}

/*
inline double Model_NC_QED::Aqed(double t) {
  return (*APHYTOOLS::aqed)(t);
}
inline double Model_NC_QED::Aqed() {
  return APHYTOOLS::aqed->Aqed(sqr(rpa.gen.Ecms()));
}
*/













