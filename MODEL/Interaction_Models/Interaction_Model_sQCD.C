#include "MODEL/Interaction_Models/Interaction_Model_sQCD.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include <stdio.h>

using namespace MODEL;
using namespace ATOOLS;
using namespace std;

Interaction_Model_sQCD::Interaction_Model_sQCD(MODEL::Model_Base * _model,
					       std::string _cplscheme,std::string _yukscheme) :
  Interaction_Model_Base("",_model,_cplscheme,_yukscheme)
{ 
  double scale = rpa.gen.CplScale();
  g1       = Kabbala(string("g_1"),
		     sqrt(4.*M_PI*ScalarFunction(string("alpha_QED"),scale)));
  g3    = Kabbala(string("g_3"),
		  sqrt(4.*M_PI*ScalarFunction(std::string("alpha_S"),scale)));
  PL    = Kabbala(string("P_L"),1.);
  PR    = Kabbala(string("P_R"),1.);
  M_I   = Kabbala(string("i"),Complex(0.,1.)); 
  root2 = Kabbala(string("\\sqrt{2}"),sqrt(2.));
  num_2    = Kabbala(string("2"),2.);    
}

void Interaction_Model_sQCD::c_FFS(std::vector<Single_Vertex>& vertex,int& vanz)
{
  Kabbala kcpl0,kcpl1;

  //quark - squark - gluino

  Flavour flgluino = Flavour(kf_Gluino);
  if (flgluino.IsOn()) {
  //uptype - sup - gluino
    for (short int i=2;i<7;i+=2) {
      Flavour flav1 = Flavour((kf_code)(i));
      int fl; 
      for (short int j=1;j<7;j++) {
	if (j<4) fl = 1000000 + 2*j;
	else     fl = 2000000 + 2*j - 6;
	Flavour flav2 = Flavour((kf_code)(fl));
	if (flav1.IsOn() && flav2.IsOn() && gen_sUp(flav2)==((i-2)/2)) {
	  vertex[vanz].in[0] = flav1;
	  vertex[vanz].in[2] = flav2;
	  vertex[vanz].in[1] = flgluino;
	  
	  kcpl0 = M_I*g3*root2*K_Z_U((i-2)/2+3,j-1);
	  kcpl1 = -M_I*g3*root2*K_Z_U((i-2)/2,j-1);
	  
	  vertex[vanz].cpl[0]  = kcpl0;
	  vertex[vanz].cpl[1]  = kcpl1;
	  vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	  
	  
	  vertex[vanz].Color.push_back(Color_Function(cf::T));     
	  
	  vertex[vanz].Color.back().SetParticleArg(1,2,0);     
	  vertex[vanz].Color.back().SetStringArg('1','2','0');     
	  
	  
	  vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("FFS",LF_Key()));
	  
	  vertex[vanz].on      = 1;
	  vertex.push_back(Single_Vertex());vanz++;
	}
      }
    }
    //downtype - sdown - gluino
    for (short int i=1;i<6;i+=2) {
      Flavour flav1 = Flavour((kf_code)(i));
      int fl;
      for (short int j=1;j<7;j++) {
	if (j<4) fl = 1000000 + 2*j -1;
	else     fl = 2000000 + 2*j -7;
	Flavour flav2 = Flavour((kf_code)(fl));
	if (flav1.IsOn() && flav2.IsOn() && gen_sDown(flav2)==((i-1)/2)) {

	  vertex[vanz].in[0] = flav1;
	  vertex[vanz].in[2] = flav2;
	  vertex[vanz].in[1] = flgluino;
	  
	  kcpl0 = M_I*g3*root2*K_Z_D((i-1)/2+3,j-1);
	  kcpl1 = -M_I*g3*root2*K_Z_D((i-1)/2,j-1);

	  
	  vertex[vanz].Color.push_back(Color_Function(cf::T));     
	  vertex[vanz].Color.back().SetParticleArg(1,2,0);     
	  vertex[vanz].Color.back().SetStringArg('1','2','0');     
	  
	  vertex[vanz].cpl[0]  = kcpl0;
	  vertex[vanz].cpl[1]  = kcpl1;
	  vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	  
	  
	  vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("FFS",LF_Key()));
	  
	  vertex[vanz].on      = 1;
	  vertex.push_back(Single_Vertex());vanz++;
	}  
      }
    }
  }
}

void Interaction_Model_sQCD::c_FFV(std::vector<Single_Vertex>& vertex,int& vanz)
{
  Kabbala kcpl0,kcpl1;
  
  //gluino - gluon - gluino
  Flavour flgluino = Flavour(kf_Gluino);
  if (flgluino.IsOn()) {      
    Flavour flgluon = Flavour(kf_gluon);
    if (flgluon.IsOn()) {
      vertex[vanz].in[0] = flgluino;
      vertex[vanz].in[1] = flgluon;
      vertex[vanz].in[2] = flgluino;
      
      kcpl0 = -g3; 
      kcpl1 = -g3;

      vertex[vanz].cpl[0]  = kcpl0;
      vertex[vanz].cpl[1]  = kcpl1;
      vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
      
      vertex[vanz].Color.push_back(Color_Function(cf::F));     
      vertex[vanz].Color.back().SetParticleArg(0,1,2);     
      vertex[vanz].Color.back().SetStringArg('0','1','2');     
      
      vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("Gamma",LF_Key()));
      vertex[vanz].Lorentz.back()->SetParticleArg(1);     
      
      vertex[vanz].on      = 1;
      vertex.push_back(Single_Vertex());vanz++;
    }   
  }
}

void Interaction_Model_sQCD::c_SSV(std::vector<Single_Vertex>& vertex,int& vanz)
{
  //sQuark - Gluon - sQuark
  
  Kabbala kcpl0,kcpl1;

  Flavour flgl = Flavour(kf_gluon); 
  if (flgl.IsOn()) {    
    kcpl0 = -g3*M_I;
    kcpl1 = kcpl0;
    for (short int l=1;l<3;l++) {
      for (short int i=1;i<7;i++) {
	int fl = l*1000000 + i;
	Flavour flav = Flavour((kf_code)(fl));
	if (flav.IsOn()) { 
	  vertex[vanz].in[0] = flav;
	  vertex[vanz].in[1] = flgl;
	  vertex[vanz].in[2] = flav;
	  
	  vertex[vanz].cpl[0]  = kcpl0;
	  vertex[vanz].cpl[1]  = kcpl1;
	  vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	  
	  
	  vertex[vanz].Color.push_back(Color_Function(cf::T));     
	  vertex[vanz].Color.back().SetParticleArg(1,2,0);     
	  vertex[vanz].Color.back().SetStringArg('1','2','0');     
	  
	  
	  vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("SSV",LF_Key()));
	  vertex[vanz].Lorentz.back()->SetParticleArg(0,2,1);     
	  
	  vertex[vanz].on      = 1;
	  vertex.push_back(Single_Vertex());vanz++;
	} 
      }
    }
  }
}

void Interaction_Model_sQCD::c_SSVV(std::vector<Single_Vertex>& vertex,int& vanz)
{
  Kabbala kcpl0,kcpl1;
  Flavour flgl = Flavour(kf_gluon); 
  
  //sQuark - Gluon - Gluon - sQuark
  if (flgl.IsOn()) {    
    for (short int l=1;l<3;l++) {
      for (short int i=1;i<7;i++) {
	int fl = l*1000000 + i;
	Flavour flav = Flavour((kf_code)(fl));
	if (flav.IsOn()) {
	  vertex[vanz].in[0] = flgl;
	  vertex[vanz].in[1] = flav.Bar();
	  vertex[vanz].in[2] = flav;
	  vertex[vanz].in[3] = flgl;
	  
	  vertex[vanz].nleg     = 4;
	  
	  kcpl0 = M_I*g3*g3;
	  kcpl1 = kcpl0;
	  
	  vertex[vanz].cpl[0]  = kcpl0;
	  vertex[vanz].cpl[1]  = kcpl1;
	  vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	  
	  
	  
	  vertex[vanz].Color.resize(2); 
	  vertex[vanz].Lorentz.resize(2); 
	  
	  vertex[vanz].Color[0]        = Color_Function(cf::T,0,2,4,'0','2','4',
							new Color_Function(cf::T,3,4,1,'3','4','1'));
	  
	  vertex[vanz].Color[1]        = Color_Function(cf::T,3,2,4,'3','2','4',
							new Color_Function(cf::T,0,4,1,'0','4','1'));
	  
	  
	  vertex[vanz].Lorentz[0]=LF_Getter::GetObject("VVSS",LF_Key());     
	  vertex[vanz].Lorentz[0]->SetParticleArg(0,3);     
	  vertex[vanz].Lorentz[1]=LF_Getter::GetObject("VVSS",LF_Key());     
	  vertex[vanz].Lorentz[1]->SetParticleArg(0,3);     
	  
	  vertex[vanz].on      = 1;
	  vertex.push_back(Single_Vertex());vanz++;
	} 
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
		 ComplexMatrixElement(std::string("Z_d"),i,j));
}  

Kabbala Interaction_Model_sQCD::K_Z_U(short int i,short int j)       
{   
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("Z^{")+string(hi)+string(hj)+string("}_U"),
		 ComplexMatrixElement(std::string("Z_u"),i,j));
}  

int Interaction_Model_sQCD::gen_sUp(Flavour fl)
{
  int gen_sUp;

  if (fl.Kfcode() == 1000002 || fl.Kfcode() == 2000002)
    gen_sUp = 0;
  if (fl.Kfcode() == 1000004 || fl.Kfcode() == 2000004)
    gen_sUp = 1;
  if (fl.Kfcode() == 1000006 || fl.Kfcode() == 2000006)
    gen_sUp = 2;

  return gen_sUp;
}

int Interaction_Model_sQCD::gen_sDown(Flavour fl)
{
  int gen_sDown;

  if (fl.Kfcode() == 1000001 || fl.Kfcode() == 2000001)
    gen_sDown = 0;
  if (fl.Kfcode() == 1000003 || fl.Kfcode() == 2000003)
    gen_sDown = 1;
  if (fl.Kfcode() == 1000005 || fl.Kfcode() == 2000005)
    gen_sDown = 2;

  return gen_sDown;
}

