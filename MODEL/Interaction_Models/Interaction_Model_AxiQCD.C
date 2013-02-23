#include "MODEL/Interaction_Models/Interaction_Model_AxiQCD.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Run_Parameter.H"

using namespace MODEL;
using namespace ATOOLS;
using namespace std;

DECLARE_GETTER(Interaction_Model_AxiQCD,"pure_AxiQCD",
	       Interaction_Model_Base,Interaction_Model_Arguments);

Interaction_Model_Base *ATOOLS::Getter
<Interaction_Model_Base,Interaction_Model_Arguments,Interaction_Model_AxiQCD>::
operator()(const Interaction_Model_Arguments &args) const
{
  return new Interaction_Model_AxiQCD
    (args.p_model,args.m_cplscheme,args.m_yukscheme);
}

void ATOOLS::Getter<Interaction_Model_Base,Interaction_Model_Arguments,
		    Interaction_Model_AxiQCD>::PrintInfo
(std::ostream &str,const size_t width) const
{ 
  str<<"The SM QCD+AxiQCD interactions only"; 
}

Interaction_Model_AxiQCD::Interaction_Model_AxiQCD(MODEL::Model_Base * _model,
					     std::string _cplscheme,std::string _yukscheme) :
  Interaction_Model_Base("pure_AxiQCD",_model,_cplscheme,_yukscheme)
{ 
  p_moqcd  = new Interaction_Model_QCD(p_model,_cplscheme,_yukscheme); 

  g3  = Kabbala(string("g_3"),sqrt(4.*M_PI*ScalarFunction(std::string("alpha_S"),rpa->gen.CplScale())));
  PL  = Kabbala(string("P_L"),1.);
  PR  = Kabbala(string("P_R"),1.);
  M_I = Kabbala(string("i"),Complex(0.,1.)); 
}

void Interaction_Model_AxiQCD::c_FFV(std::vector<Single_Vertex>& vertex,int & vanz)
{
  p_moqcd->c_FFV(vertex,vanz);

  Kabbala kcpl0 = -g3*M_I;
  Kabbala kcpl1 = kcpl0;

  for (short int i=1;i<=6;i++) {
    Flavour flav = Flavour((kf_code)(i));
    if (flav.Strong() && flav.IsOn() && Flavour(61).IsOn()) { 
      vertex[vanz].in[0]         = flav;
      vertex[vanz].in[1]         = Flavour(61);
      vertex[vanz].in[2]         = flav;

      vertex[vanz].cpl[0]        = kcpl0;
      vertex[vanz].cpl[1]        = kcpl1;
      vertex[vanz].Str           = (kcpl0*PR-kcpl1*PL).String();      

      
      vertex[vanz].Color.push_back(Color_Function(cf::T,1,2,0,'1','2','0'));     

      vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("Gamma",LF_Key()));     
      vertex[vanz].Lorentz.back()->SetParticleArg(1);     
                  
      vertex[vanz].on            = 1;
      vertex.push_back(Single_Vertex());vanz++;
    } 
  }
}

void Interaction_Model_AxiQCD::c_VVV(std::vector<Single_Vertex>& vertex,int& vanz)
{
  p_moqcd->c_VVV(vertex,vanz);

  Kabbala kcpl0 = -g3;
  Kabbala kcpl1 = kcpl0; 
  
  if (Flavour(kf_gluon).IsOn() && Flavour(61).IsOn()) {

    vertex[vanz].in[0] = Flavour(kf_gluon);
    vertex[vanz].in[1] = Flavour(61);
    vertex[vanz].in[2] = Flavour(61);

    vertex[vanz].cpl[0]        = kcpl0;
    vertex[vanz].cpl[1]        = kcpl1;
    vertex[vanz].Str           = (kcpl0*PR+kcpl1*PL).String();
    
    
    vertex[vanz].Color.push_back(Color_Function(cf::F));     
    vertex[vanz].Color.back().SetParticleArg(0,2,1);     
    vertex[vanz].Color.back().SetStringArg('0','2','1');     
    
    vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("Gauge3",LF_Key()));     
    vertex[vanz].Lorentz.back()->SetParticleArg(0,1,2);     
    
    vertex[vanz].on            = 1;
    vertex.push_back(Single_Vertex());vanz++;
  }
}

void Interaction_Model_AxiQCD::c_VVVV(std::vector<Single_Vertex>& vertex,int& vanz)
{
  p_moqcd->c_VVVV(vertex,vanz);

  Kabbala kcpl0 = -M_I*g3*g3; 
  Kabbala kcpl1 = kcpl0; 
  
  if (Flavour(kf_gluon).IsOn() && Flavour(61).IsOn()) {

    vertex[vanz].in[0] = Flavour(kf_gluon);
    vertex[vanz].in[1] = Flavour(kf_gluon);
    vertex[vanz].in[2] = Flavour(61);
    vertex[vanz].in[3] = Flavour(61);

    vertex[vanz].nleg            = 4;
    vertex[vanz].cpl[0]          = kcpl0;
    vertex[vanz].cpl[1]          = kcpl1;
    vertex[vanz].Str             = (kcpl0*PR+kcpl1*PL).String();
    
    vertex[vanz].Color.resize(3);
    vertex[vanz].Lorentz.resize(3); 
    
    vertex[vanz].Color[0]        = Color_Function(cf::F,0,2,4,'0','2','4',
						  new Color_Function(cf::F,1,3,4,'1','3','4'));
    vertex[vanz].Lorentz[0]=LF_Getter::GetObject("Gluon4",LF_Key());
    vertex[vanz].Lorentz[0]->SetParticleArg(0,1,2,3);     
    
    vertex[vanz].Color[1]        = Color_Function(cf::F,0,3,4,'0','3','4',
						  new Color_Function(cf::F,1,2,4,'1','2','4'));
    vertex[vanz].Lorentz[1]=LF_Getter::GetObject("Gluon4",LF_Key());
    vertex[vanz].Lorentz[1]->SetParticleArg(0,1,3,2);     
    
    vertex[vanz].Color[2]        = Color_Function(cf::F,0,1,4,'0','1','4',
						  new Color_Function(cf::F,3,2,4,'3','2','4')); 
    vertex[vanz].Lorentz[2]=LF_Getter::GetObject("Gluon4",LF_Key());     
    vertex[vanz].Lorentz[2]->SetParticleArg(0,3,1,2);     
    
    vertex[vanz].on              = 1;
    vertex.push_back(Single_Vertex());vanz++;




    vertex[vanz].in[0] = Flavour(61);
    vertex[vanz].in[1] = Flavour(61);
    vertex[vanz].in[2] = Flavour(61);
    vertex[vanz].in[3] = Flavour(61);

    vertex[vanz].nleg            = 4;
    vertex[vanz].cpl[0]          = kcpl0;
    vertex[vanz].cpl[1]          = kcpl1;
    vertex[vanz].Str             = (kcpl0*PR+kcpl1*PL).String();
    
    vertex[vanz].Color.resize(3);
    vertex[vanz].Lorentz.resize(3); 
    
    vertex[vanz].Color[0]        = Color_Function(cf::F,0,2,4,'0','2','4',
						  new Color_Function(cf::F,1,3,4,'1','3','4'));
    vertex[vanz].Lorentz[0]=LF_Getter::GetObject("Gluon4",LF_Key());
    vertex[vanz].Lorentz[0]->SetParticleArg(0,1,2,3);     
    
    vertex[vanz].Color[1]        = Color_Function(cf::F,0,3,4,'0','3','4',
						  new Color_Function(cf::F,1,2,4,'1','2','4'));
    vertex[vanz].Lorentz[1]=LF_Getter::GetObject("Gluon4",LF_Key());
    vertex[vanz].Lorentz[1]->SetParticleArg(0,1,3,2);     
    
    vertex[vanz].Color[2]        = Color_Function(cf::F,0,1,4,'0','1','4',
						  new Color_Function(cf::F,3,2,4,'3','2','4')); 
    vertex[vanz].Lorentz[2]=LF_Getter::GetObject("Gluon4",LF_Key());     
    vertex[vanz].Lorentz[2]->SetParticleArg(0,3,1,2);     
    
    vertex[vanz].on              = 1;
    vertex.push_back(Single_Vertex());vanz++;
  }  
}

