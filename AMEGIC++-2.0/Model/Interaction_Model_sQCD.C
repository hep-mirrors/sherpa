#include "Interaction_Model_sQCD.H"
#include "MathTools.H"
#include "Message.H"
#include <stdio.h>

using namespace AMEGIC;
using namespace ATOOLS;
using namespace std;

Interaction_Model_sQCD::Interaction_Model_sQCD(MODEL::Model_Base * _model,
					       std::string _cplscheme,std::string _yukscheme) :
  Interaction_Model_Base(_model,_cplscheme,_yukscheme)
{ 
  double Ecms2 = sqr(rpa.gen.Ecms());
  g1       = Kabbala(string("g_1"),
		     sqrt(4.*M_PI*ScalarFunction(string("alpha_QED"),Ecms2)));
  g3    = Kabbala(string("g_3"),
		  sqrt(4.*M_PI*ScalarFunction(std::string("alpha_S"),Ecms2)));
  PL    = Kabbala(string("P_L"),1.);
  PR    = Kabbala(string("P_R"),1.);
  M_I   = Kabbala(string("i"),Complex(0.,1.)); 
  root2 = Kabbala(string("\\sqrt{2}"),sqrt(2.));
  num_2    = Kabbala(string("2"),2.);    
}

void Interaction_Model_sQCD::c_FFS(Single_Vertex* vertex,int& vanz)
{
  Kabbala kcpl0,kcpl1;

  //quark - squark - gluino

  Flavour flgluino = Flavour(kf::code(47));
  if (flgluino.IsOn()) {
  //uptype - sup - gluino
    for (short int i=2;i<7;i+=2) {
      Flavour flav1 = Flavour(kf::code(i));
      for (short int j=51;j<57;j++) {
	Flavour flav2 = Flavour(kf::code(j));
	if (flav1.IsOn() && flav2.IsOn()) {
	vertex[vanz].in[0] = flav1;
	vertex[vanz].in[1] = flav2;
	vertex[vanz].in[2] = flgluino;
	
	kcpl0 = M_I*g3*root2*K_Z_U((i-2)/2+3,j-51);
	kcpl1 = -M_I*g3*root2*K_Z_U((i-2)/2,j-51);
	
	vertex[vanz].cpl[0]  = kcpl0.Value();
	vertex[vanz].cpl[1]  = kcpl1.Value();
	vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;

	vertex[vanz].ncf   = 1;
	vertex[vanz].Color = new Color_Function(cf::T);     
	vertex[vanz].Color->SetParticleArg(2,1,0);     
	vertex[vanz].Color->SetStringArg('2','1','0');     
	
	vertex[vanz].nlf     = 1;
	vertex[vanz].Lorentz = new Lorentz_Function(lf::FFS);
	
	vertex[vanz].on      = 1;
	vanz++;
      }  
      }
    }
  //downtype - sdown - gluino
    for (short int i=1;i<6;i+=2) {
      Flavour flav1 = Flavour(kf::code(i));
      for (short int j=61;j<67;j++) {
	Flavour flav2 = Flavour(kf::code(j));
	if (flav1.IsOn() && flav2.IsOn()) {
	vertex[vanz].in[0] = flav1;
	vertex[vanz].in[1] = flav2;
	vertex[vanz].in[2] = flgluino;
	
	kcpl0 = M_I*g3*root2*K_Z_D((i-1)/2+3,j-61);
	kcpl1 = -M_I*g3*root2*K_Z_D((i-1)/2,j-61);
	
	vertex[vanz].ncf   = 1;
	vertex[vanz].Color = new Color_Function(cf::T);     
	vertex[vanz].Color->SetParticleArg(2,1,0);     
	vertex[vanz].Color->SetStringArg('2','1','0');     

	vertex[vanz].cpl[0]  = kcpl0.Value();
	vertex[vanz].cpl[1]  = kcpl1.Value();
	vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;
	
	vertex[vanz].nlf     = 1;
	vertex[vanz].Lorentz = new Lorentz_Function(lf::FFS);
	
	vertex[vanz].on      = 1;
	vanz++;
	}  
      }
    }
  }
}

void Interaction_Model_sQCD::c_FFV(Single_Vertex* vertex,int& vanz)
{
  Kabbala kcpl0,kcpl1;
  
  //gluino - gluon - gluino
  Flavour flgluino = Flavour(kf::code(47));
  if (flgluino.IsOn()) {      
    Flavour flgluon = Flavour(kf::gluon);
    if (flgluon.IsOn()) {
      vertex[vanz].in[0] = flgluino;
      vertex[vanz].in[1] = flgluon;
      vertex[vanz].in[2] = flgluino;
      
      kcpl0 = -g3; 
      kcpl1 = -g3;
      
      vertex[vanz].cpl[0]  = kcpl0.Value();
      vertex[vanz].cpl[1]  = kcpl1.Value();
      vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
      vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;
      
      vertex[vanz].ncf   = 1;
      vertex[vanz].Color = new Color_Function(cf::F);     
      vertex[vanz].Color->SetParticleArg(0,1,2);     
      vertex[vanz].Color->SetStringArg('0','1','2');     
      
      vertex[vanz].nlf     = 1;
      vertex[vanz].Lorentz = new Lorentz_Function(lf::Gamma);
      vertex[vanz].Lorentz->SetParticleArg(1);     
      
      vertex[vanz].on      = 1;
      vanz++;
    }   
  }
}

void Interaction_Model_sQCD::c_SSV(Single_Vertex* vertex,int& vanz)
{
  //sQuark - Gluon - sQuark
  
  Kabbala kcpl0,kcpl1;

  Flavour flgl = Flavour(kf::gluon); 
  if (flgl.IsOn()) {    
    //uptypes 
    kcpl0 = g3*M_I;
    kcpl1 = kcpl0;
    for (short int i=51;i<57;i++) {
      Flavour flav = Flavour(kf::code(i));
      if (flav.IsOn()) { 
	vertex[vanz].in[0] = flav;
	vertex[vanz].in[1] = flgl;
	vertex[vanz].in[2] = flav;	
	vertex[vanz].cpl[0]  = kcpl0.Value();
	vertex[vanz].cpl[1]  = kcpl1.Value();
	vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;
	
	vertex[vanz].ncf   = 1;
	vertex[vanz].Color = new Color_Function(cf::T);     
	vertex[vanz].Color->SetParticleArg(1,0,2);     
	vertex[vanz].Color->SetStringArg('1','0','2');     
	
	vertex[vanz].nlf     = 1;
	vertex[vanz].Lorentz = new Lorentz_Function(lf::SSV);
	vertex[vanz].Lorentz->SetParticleArg(0,2,1);     
	
	vertex[vanz].on      = 1;
	vanz++;
      } 
    }
    //downtypes 
    kcpl0 = -g3*M_I;
    kcpl1 = kcpl0;
    for (short int i=61;i<67;i++) {
      Flavour flav = Flavour(kf::code(i));
      if (flav.IsOn()) { 
	vertex[vanz].in[0] = flav;
	vertex[vanz].in[1] = flgl;
	vertex[vanz].in[2] = flav;
	
	vertex[vanz].cpl[0]  = kcpl0.Value();
	vertex[vanz].cpl[1]  = kcpl1.Value();
	vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;
	
	vertex[vanz].ncf   = 1;
	vertex[vanz].Color = new Color_Function(cf::T);     
	vertex[vanz].Color->SetParticleArg(1,0,2);     
	vertex[vanz].Color->SetStringArg('1','0','2');     
	
	vertex[vanz].nlf     = 1;
	vertex[vanz].Lorentz = new Lorentz_Function(lf::SSV);
	vertex[vanz].Lorentz->SetParticleArg(0,2,1);     
	
	vertex[vanz].on      = 1;
	vanz++;
      } 
    }
  }
}

void Interaction_Model_sQCD::c_SSVV(Single_Vertex* vertex,int& vanz)
{
  Kabbala kcpl0,kcpl1;
  Flavour flgl = Flavour(kf::gluon); 
  Flavour flph = Flavour(kf::photon); 
  
  //sQuark - Gluon - Gluon - sQuark
  if (flgl.IsOn()) {    
    for (short int i=51;i<67;i++) {
      if (i==57) i=61;
      Flavour flav = Flavour(kf::code(i));
      if (flav.IsOn()) {
	vertex[vanz].in[0] = flgl;
	vertex[vanz].in[1] = flav.Bar();
	vertex[vanz].in[2] = flav;
	vertex[vanz].in[3] = flgl;
	
	vertex[vanz].nleg     = 4;
	
	kcpl0 = M_I*g3*g3;
	kcpl1 = kcpl0;
		
	vertex[vanz].cpl[0]  = kcpl0.Value();
	vertex[vanz].cpl[1]  = kcpl1.Value();
	vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;
	
	vertex[vanz].ncf     = 2;
	vertex[vanz].nlf     = 2;
	vertex[vanz].Color   = new Color_Function[2]; 
	vertex[vanz].Lorentz = new Lorentz_Function[2]; 

	vertex[vanz].Color[0]        = Color_Function(cf::T,0,1,4,'0','1','4',
						      new Color_Function(cf::T,3,4,2,'3','4','2');
		
	vertex[vanz].Color[1]        = Color_Function(cf::T,3,1,4,'3','1','4',
						      new Color_Function(cf::T,0,4,2,'0','4','2');
						      
						      
	vertex[vanz].Lorentz[0] = Lorentz_Function(lf::VVSS);     
	vertex[vanz].Lorentz[0].SetParticleArg(0,3);     
	vertex[vanz].Lorentz[1] = Lorentz_Function(lf::VVSS);     
	vertex[vanz].Lorentz[1].SetParticleArg(0,3);     
	
	vertex[vanz].on      = 1;
	vanz++;
      } 
    }
  }
 
  //sQuark - Gluon - Photon - sQuark
  if (flgl.IsOn() && flph.IsOn()) {    
    for (short int i=51;i<67;i++) {
      if (i==57) i=61;
      Flavour flav = Flavour(kf::code(i));
      if (flav.IsOn()) {
	vertex[vanz].in[0] = flgl;
	vertex[vanz].in[1] = flav.Bar();
	vertex[vanz].in[2] = flav;
	vertex[vanz].in[3] = flph;
	
	vertex[vanz].nleg     = 4;
	
	Kabbala charge = Kabbala(string("Q_{"+flav.TexName()+"}"),flav.Charge());
	
	kcpl0 = M_I*num_2*charge*g1*g3;
	kcpl1 = kcpl0;
	
	vertex[vanz].cpl[0]  = kcpl0.Value();
	vertex[vanz].cpl[1]  = kcpl1.Value();
	vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;
	
	vertex[vanz].ncf     = 1;
	vertex[vanz].Color   = new Color_Function(cf::T); 
	vertex[vanz].Color->SetParticleArg(0,1,2);     
	vertex[vanz].Color->SetStringArg('0','1','2');     
	
	vertex[vanz].nlf     = 1;
	vertex[vanz].Lorentz = new Lorentz_Function(lf::VVSS);     
	vertex[vanz].Lorentz->SetParticleArg(0,3);     
	
	vertex[vanz].on      = 1;
	vanz++;
      } 
    }
  }
}

Kabbala Interaction_Model_sQCD::K_Z_D(short int i,short int j)       
{   
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("Z^{")+string(hi)+string(hj)+string("}_D"),
		 p_model->ComplexMatrixElement(std::string("Z_d"),i,j));
}  

Kabbala Interaction_Model_sQCD::K_Z_U(short int i,short int j)       
{   
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("Z^{")+string(hi)+string(hj)+string("}_U"),
		 p_model->ComplexMatrixElement(std::string("Z_u"),i,j));
}  

