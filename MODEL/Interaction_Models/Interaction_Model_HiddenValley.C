#include "MODEL/Interaction_Models/Interaction_Model_HiddenValley.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include <stdio.h>


using namespace MODEL;
using namespace ATOOLS;
using namespace std;

DECLARE_GETTER(Interaction_Model_HiddenValley_Getter,"SM+HiddenValley",
	       Interaction_Model_Base,Interaction_Model_Arguments);

Interaction_Model_Base *Interaction_Model_HiddenValley_Getter::
operator()(const Interaction_Model_Arguments &args) const
{
  return new Interaction_Model_HiddenValley
    (args.p_model,args.m_cplscheme,args.m_yukscheme);
}

void Interaction_Model_HiddenValley_Getter::PrintInfo
(std::ostream &str,const size_t width) const
{ 
  str<<"The Standard Model plus hidden sector"; 
}

Interaction_Model_HiddenValley::Interaction_Model_HiddenValley(MODEL::Model_Base * _model,
					   std::string _cplscheme,std::string _yukscheme) :
  Interaction_Model_Base("HiddenValley",_model,_cplscheme,_yukscheme)
{ 
  p_mosm  = new Interaction_Model_SM(p_model,_cplscheme,_yukscheme); 

  g1       = Kabbala(string("g_1"),
		     sqrt(4.*M_PI*ScalarFunction(std::string("alpha_QED"),rpa.gen.CplScale())));
  g2       = Kabbala(string("g_1/\\sin\\theta_W"), 
		     g1.Value()/sqrt(ScalarConstant(std::string("sin2_thetaW"))));
  gD3  = Kabbala(string("gD_3"),sqrt(4.*M_PI*ScalarFunction(std::string("alpha_HV"),rpa.gen.CplScale())));
  
  sintW    = Kabbala(std::string("\\sin\\theta_W"),
		     sqrt(ScalarConstant(std::string("sin2_thetaW"))));
  costW    = Kabbala(std::string("\\cos\\theta_W"),
		     sqrt(1.-ScalarConstant(std::string("sin2_thetaW"))));

  PL       = Kabbala(string("P_L"),1.);
  PR       = Kabbala(string("P_R"),1.);
  M_I      = Kabbala(string("i"),Complex(0.,1.));
  sqrt_HV_Nc    = Kabbala(string("\\sqrt{HV_Nc}"),sqrt(ScalarConstant(std::string("HV_NC"))));
  sqrt_3        = Kabbala(string("\\sqrt{3}"),sqrt(3.));
  
}

void Interaction_Model_HiddenValley::c_FFV(std::vector<Single_Vertex>& vertex,int& vanz)
{
  p_mosm->c_FFV(vertex,vanz);

  Kabbala kcpl0,kcpl1;
  Flavour flZp(9900023);
  //visible sector
  for (short int i=1;i<17;i++) {
    if (i==7) i=11;
    Flavour SMferm    = Flavour((kf_code)(i));
    Kabbala Q   = Kabbala(string("Q_{")+SMferm.TexName()+string("}"),SMferm.Charge());
    Kabbala TR  = Kabbala(string("T_{")+SMferm.TexName()+string("}"),SMferm.IsoWeak());
    
    kcpl0             = M_I/costW*Q*sintW*sintW*g2;
    kcpl1             = -M_I/costW*(TR-Q*sintW*sintW)*g2;
    
    vertex[vanz].in[0]     = SMferm;
    vertex[vanz].in[1]     = flZp;
    vertex[vanz].in[2]     = SMferm;
    vertex[vanz].cpl[0]    = kcpl0;
    vertex[vanz].cpl[1]    = kcpl1;
    vertex[vanz].Str       = (kcpl0*PR+kcpl1*PL).String();
    
    
    if (SMferm.Strong()) {
      vertex[vanz].Color.push_back(Color_Function(cf::D));;     
      vertex[vanz].Color.back().SetParticleArg(0,2);     
      vertex[vanz].Color.back().SetStringArg('0','2');     
    }
    else 
      vertex[vanz].Color.push_back(Color_Function(cf::None));;
    
    
    vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("Gamma",LF_Key()));
    vertex[vanz].Lorentz.back()->SetParticleArg(1);     
    
    vertex[vanz].on     = 1;
    vertex.push_back(Single_Vertex());vanz++;
  }
  
  //hidden sector
  //messenger Zp
  for (short int i=1;i<3;i++) {
    Flavour HVferm    = Flavour((kf_code)(9900000+i));
    Kabbala Q   = Kabbala(string("Q_{")+HVferm.TexName()+string("}"),HVferm.Charge());
    Kabbala TR  = Kabbala(string("T_{")+HVferm.TexName()+string("}"),HVferm.IsoWeak());
    
    Kabbala col_fac = sqrt_HV_Nc/sqrt_3;

    kcpl0 = M_I*col_fac/costW*Q*sintW*sintW*g2;
    kcpl1 = -M_I*col_fac/costW*(TR-Q*sintW*sintW)*g2;
    
    vertex[vanz].in[0]     = HVferm;
    vertex[vanz].in[1]     = flZp;
    vertex[vanz].in[2]     = HVferm;
    vertex[vanz].cpl[0]    = kcpl0;
    vertex[vanz].cpl[1]    = kcpl1;
    vertex[vanz].Str       = (kcpl0*PR+kcpl1*PL).String();
    
    //might be changed to none
    //adjust coupling instead
    if (HVferm.Strong()) {
      vertex[vanz].Color.push_back(Color_Function(cf::D));;     
      vertex[vanz].Color.back().SetParticleArg(0,2);     
      vertex[vanz].Color.back().SetStringArg('0','2');     
    }
    else 
      vertex[vanz].Color.push_back(Color_Function(cf::None));;
    
    
    vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("Gamma",LF_Key()));
    vertex[vanz].Lorentz.back()->SetParticleArg(1);     
    
    vertex[vanz].on     = 1;
    vertex.push_back(Single_Vertex());vanz++;
  }
  //hidden gluon
  Flavour flgD(9900021);
  if (flgD.IsOn()) {
    for (short int i=1;i<3;i++) {
      Flavour HVferm    = Flavour((kf_code)(9900000+i));
      
      kcpl0 = kcpl1 = M_I*gD3;
      
      vertex[vanz].in[0]     = HVferm;
      vertex[vanz].in[1]     = flgD;
      vertex[vanz].in[2]     = HVferm;
      vertex[vanz].cpl[0]    = kcpl0;
      vertex[vanz].cpl[1]    = kcpl1;
      vertex[vanz].Str       = (kcpl0*PR+kcpl1*PL).String();
      
      //might be changed to none
      //adjust coupling instead
      if (HVferm.Strong()) {
	vertex[vanz].Color.push_back(Color_Function(cf::D));;     
	vertex[vanz].Color.back().SetParticleArg(0,2);     
	vertex[vanz].Color.back().SetStringArg('0','2');     
      }
      else 
	vertex[vanz].Color.push_back(Color_Function(cf::None));;
      
      
      vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("Gamma",LF_Key()));
      vertex[vanz].Lorentz.back()->SetParticleArg(1);     
      
      vertex[vanz].on     = 1;
      vertex.push_back(Single_Vertex());vanz++;
    }
  }
}

void Interaction_Model_HiddenValley::c_VVV(std::vector<Single_Vertex>& vertex,int& vanz)
{
  Kabbala kcpl0 = -gD3;
  Kabbala kcpl1 = kcpl0; 
  
  Flavour flgD = Flavour(9900021);
  if (flgD.IsOn()) {

    for (short int i=0;i<3;i++) vertex[vanz].in[i] = flgD;

  vertex[vanz].cpl[0]        = kcpl0;
  vertex[vanz].cpl[1]        = kcpl1;
  vertex[vanz].Str           = (kcpl0*PR+kcpl1*PL).String();

  
  vertex[vanz].Color.push_back(Color_Function(cf::F));;     
  vertex[vanz].Color.back().SetParticleArg(0,2,1);     
  vertex[vanz].Color.back().SetStringArg('0','2','1');     

  vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("Gauge3",LF_Key()));     
  vertex[vanz].Lorentz.back()->SetParticleArg(0,1,2);     

  vertex[vanz].on            = 1;
  vertex.push_back(Single_Vertex());vanz++;
  }
  //SM triple vectors
  p_mosm->c_VVV(vertex,vanz);
}
void Interaction_Model_HiddenValley::c_VVVV(std::vector<Single_Vertex>& vertex,int& vanz)
{
  p_mosm->c_VVVV(vertex,vanz);
}

void Interaction_Model_HiddenValley::c_FFS(std::vector<Single_Vertex>& vertex,int& vanz)  
{ 
  p_mosm->c_FFS(vertex,vanz); 
}
void Interaction_Model_HiddenValley::c_VVS(std::vector<Single_Vertex>& vertex,int& vanz)  
{ 
  p_mosm->c_VVS(vertex,vanz); 
}
void Interaction_Model_HiddenValley::c_SSS(std::vector<Single_Vertex>& vertex,int& vanz)  
{ 
  p_mosm->c_SSS(vertex,vanz); 
}
void Interaction_Model_HiddenValley::c_SSVV(std::vector<Single_Vertex>& vertex,int& vanz) 
{ 
  p_mosm->c_SSVV(vertex,vanz); 
}
void Interaction_Model_HiddenValley::c_SSSS(std::vector<Single_Vertex>& vertex,int& vanz) 
{ 
  p_mosm->c_SSSS(vertex,vanz); 
}

Interaction_Model_HiddenValley::~Interaction_Model_HiddenValley()
{
  delete p_mosm;
}
