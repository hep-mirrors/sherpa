#include "MODEL/Interaction_Models/Interaction_Model_MSSM_EHC.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include <stdio.h>


using namespace MODEL;
using namespace ATOOLS;
using namespace std;

DECLARE_GETTER(Interaction_Model_MSSM_EHC_Getter,"MSSM+EHC",
	       Interaction_Model_Base,Interaction_Model_Arguments);

Interaction_Model_Base *Interaction_Model_MSSM_EHC_Getter::
operator()(const Interaction_Model_Arguments &args) const
{
  return new Interaction_Model_MSSM_EHC
    (args.p_model,args.m_cplscheme,args.m_yukscheme);
}

void Interaction_Model_MSSM_EHC_Getter::PrintInfo
(std::ostream &str,const size_t width) const
{ 
  str<<"The MSSM + Effective Higgs Gluon Coupling"; 
}

Interaction_Model_MSSM_EHC::Interaction_Model_MSSM_EHC(MODEL::Model_Base * _model,
					       std::string _cplscheme,std::string _yukscheme) :
  Interaction_Model_Base("MSSM+EHC",_model,_cplscheme,_yukscheme)
{ 
  m_loops=true;
  p_mothdm    = new Interaction_Model_THDM(p_model,_cplscheme,_yukscheme); 
  p_moinos    = new Interaction_Model_Inos(p_model,_cplscheme,_yukscheme); 
  p_moslepton = new Interaction_Model_sLepton_EW(p_model,_cplscheme,_yukscheme); 
  p_mosqcd    = new Interaction_Model_sQCD(p_model,_cplscheme,_yukscheme); 
  p_mosquark  = new Interaction_Model_sQuark_EW(p_model,_cplscheme,_yukscheme); 
  p_moslesqu  = new Interaction_Model_sLepton_sQuark(p_model,_cplscheme,_yukscheme); 

  double scale = rpa->gen.CplScale();
  g1    = Kabbala(string("g_1"),
		  sqrt(4.*M_PI*ScalarFunction(std::string("alpha_QED"),scale)));
  g2    = Kabbala(string("g_1/\\sin\\theta_W"), 
		  g1.Value()/sqrt(ScalarConstant(std::string("sin2_thetaW"))));
  g3  = Kabbala(string("g_3"),sqrt(4.*M_PI*ScalarFunction(std::string("alpha_S"),scale)));
  sintW = Kabbala(std::string("\\sin\\theta_W"),
		  sqrt(ScalarConstant(std::string("sin2_thetaW"))));
  costW = Kabbala(std::string("\\cos\\theta_W"),
		  sqrt(1.-ScalarConstant(std::string("sin2_thetaW"))));
  PL    = Kabbala(string("P_L"),1.);
  PR    = Kabbala(string("P_R"),1.);
  M_I   = Kabbala(string("i"),Complex(0.,1.));
  root2 = Kabbala(string("\\sqrt{2}"),sqrt(2.));
  vev   = Kabbala(string("v_{EW}"),ScalarConstant(std::string("vev")));
}

void Interaction_Model_MSSM_EHC::c_FFV(std::vector<Single_Vertex>& vertex,int& vanz)
{
  p_mothdm->c_FFV(vertex,vanz);
  p_moinos->c_FFV(vertex,vanz);
  p_moslepton->c_FFV(vertex,vanz);
  p_mosqcd->c_FFV(vertex,vanz);
  p_mosquark->c_FFV(vertex,vanz);
}

void Interaction_Model_MSSM_EHC::c_VVV(std::vector<Single_Vertex>& vertex,int& vanz)
{
  p_mothdm->c_VVV(vertex,vanz);
  p_moinos->c_VVV(vertex,vanz);
  p_moslepton->c_VVV(vertex,vanz);
  p_mosqcd->c_VVV(vertex,vanz);
  p_mosquark->c_VVV(vertex,vanz);
}

void Interaction_Model_MSSM_EHC::c_VVVV(std::vector<Single_Vertex>& vertex,int& vanz)
{
  p_mothdm->c_VVVV(vertex,vanz);
  p_moinos->c_VVVV(vertex,vanz);
  p_moslepton->c_VVVV(vertex,vanz);
  p_mosqcd->c_VVVV(vertex,vanz);
  p_mosquark->c_VVVV(vertex,vanz);
  
  Kabbala kcpl0,kcpl1,ghgg; 
  Flavour flg(kf_gluon);
  
  for (int i=25;i<36;i+=10) {
    Flavour flh = Flavour((kf_code)(i));
    if(flh.IsOn() && flg.IsOn()) {
    
      // 3 gluon higgs
      vertex[vanz].nleg    = 4;
      for (short int i=0;i<3;i++)
	vertex[vanz].in[i] = flg;
      vertex[vanz].in[3] = flh;
      
      ghgg  = Kabbala(string("g_{")+flh.TexName()+string("gg}"),ScalarConstant(flh.IDName()+string("_gg_fac"))*
		      ScalarFunction(string("alpha_S"),sqr(flh.Mass()))/(2.*M_PI)/ScalarConstant(std::string("vev")));
      
      kcpl0 = g3*ghgg; 
      kcpl1 = kcpl0; 
      
      vertex[vanz].cpl[0]  = kcpl0;
      vertex[vanz].cpl[1]  = kcpl1;
      vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
      
      vertex[vanz].Color.push_back(Color_Function(cf::F));     
      vertex[vanz].Color.back().SetParticleArg(0,2,1);     
      vertex[vanz].Color.back().SetStringArg('0','2','1');     
      
      vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("Box",LF_Key()));     
      vertex[vanz].Lorentz.back()->SetParticleArg(0,1,2);     
      
      vertex[vanz].on      = 1;
      vertex.push_back(Single_Vertex());vanz++;
    }
    Flavour flsh(kf_shgluon);
    kcpl0 = M_I*g3*g3*ghgg; 
    kcpl1 = kcpl0; 
    
    if(flsh.IsOn()) {
      for (short int i=0;i<3;i++) vertex[vanz].in[i] = flg;
      vertex[vanz].in[3] = flsh;
      
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
      vertex[vanz].t               = 1;
      vertex.push_back(Single_Vertex());vanz++;
    }
  }
  
  Flavour flA = Flavour((kf_code)(kf_A0));
  if(flA.IsOn() && flg.IsOn()) {
    
    // 3 gluon higgs
    vertex[vanz].nleg    = 4;
    for (short int i=0;i<3;i++)
      vertex[vanz].in[i] = flg;
    vertex[vanz].in[3] = flA;
    
    ghgg  = Kabbala(string("g_{")+flA.TexName()+std::string("gg}"),ScalarConstant(string("A0_gg_fac"))*
		    ScalarFunction(string("alpha_S"),sqr(flA.Mass()))/(2.*M_PI)/ScalarConstant(string("vev")));
    
    kcpl0 = -g3*ghgg; 
    kcpl1 = kcpl0; 
    
    vertex[vanz].cpl[0]  = kcpl0;
    vertex[vanz].cpl[1]  = kcpl1;
    vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
    
    vertex[vanz].Color.push_back(Color_Function(cf::F));     
    vertex[vanz].Color.back().SetParticleArg(0,2,1);     
    vertex[vanz].Color.back().SetStringArg('0','2','1');     
    
    vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("PseudoBox",LF_Key()));     
    vertex[vanz].Lorentz.back()->SetParticleArg(0,1,2);     
    
    vertex[vanz].on      = 1;
    vertex.push_back(Single_Vertex());vanz++;
  }
  
}

void Interaction_Model_MSSM_EHC::c_FFS(std::vector<Single_Vertex>& vertex,int& vanz)  
{ 
  p_mothdm->c_FFS(vertex,vanz);
  p_moinos->c_FFS(vertex,vanz);
  p_moslepton->c_FFS(vertex,vanz);
  p_mosqcd->c_FFS(vertex,vanz);
  p_mosquark->c_FFS(vertex,vanz);
}

void Interaction_Model_MSSM_EHC::c_VVS(std::vector<Single_Vertex>& vertex,int& vanz)  
{ 
  p_mothdm->c_VVS(vertex,vanz);
  p_moinos->c_VVS(vertex,vanz);
  p_moslepton->c_VVS(vertex,vanz);
  p_mosqcd->c_VVS(vertex,vanz);
  p_mosquark->c_VVS(vertex,vanz);

  Kabbala kcpl0,kcpl1,ghgg;  
  Kabbala num_2 = Kabbala(string("2"),2.);  
  
  Flavour flg(kf_gluon);
  for (int i=25;i<36;i+=10) {
    Flavour flh = Flavour((kf_code)(i));
    if (flh.IsOn() && flg.IsOn()) {
      // Gluon h Gluon
      vertex[vanz].in[0] = flg;
      vertex[vanz].in[1] = flh;
      vertex[vanz].in[2] = flg;
      
      ghgg  = Kabbala(string("g_{")+flh.TexName()+std::string("gg}"),ScalarConstant(flh.IDName()+string("_gg_fac"))*
		      ScalarFunction(string("alpha_S"),sqr(flh.Mass()))/(2.*M_PI)/ScalarConstant(string("vev")));
      
      kcpl0 = M_I*ghgg;
      kcpl1 = kcpl0;
      
      vertex[vanz].cpl[0]  = kcpl0;
      vertex[vanz].cpl[1]  = kcpl0;
      vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
      
      vertex[vanz].Color.push_back(Color_Function(cf::G));     
      vertex[vanz].Color.back().SetParticleArg(0,2);     
      vertex[vanz].Color.back().SetStringArg('0','2');     
      
      vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("Triangle",LF_Key()));     
      vertex[vanz].Lorentz.back()->SetParticleArg(0,2);     
      
      vertex[vanz].on      = 1;
      vertex.push_back(Single_Vertex());vanz++;
    }
    Flavour flsh(kf_shgluon);
    // gluon h shgluon
    if (flg.IsOn() && flsh.IsOn()) {
      vertex[vanz].in[2] = flg;
      vertex[vanz].in[0] = flsh;
      vertex[vanz].in[1] = flh;
      
      kcpl0 = M_I;
      kcpl1 = kcpl0;
      
      vertex[vanz].cpl[0]  = kcpl0;
      vertex[vanz].cpl[1]  = kcpl0;
      vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
      
      vertex[vanz].Color.push_back(Color_Function(cf::G));     
      vertex[vanz].Color.back().SetParticleArg(0,2);     
      vertex[vanz].Color.back().SetStringArg('0','2');     
      
      vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("C4GS",LF_Key()));     
      vertex[vanz].Lorentz.back()->SetParticleArg(0,2);     
      
      vertex[vanz].on      = 1;
      vertex[vanz].t       = -1;
      vertex.push_back(Single_Vertex());vanz++;
    }
  } 

  Flavour flA((kf_code)(kf_A0));
  if (flA.IsOn() && flg.IsOn()) {
    
    // Gluon A0 Gluon
    vertex[vanz].in[0] = flg;
    vertex[vanz].in[1] = flA;
    vertex[vanz].in[2] = flg;
    
    ghgg  = Kabbala(std::string("g_{")+flA.TexName()+std::string("gg}"),ScalarConstant(std::string("A0_gg_fac"))*
		    ScalarFunction(std::string("alpha_S"),sqr(flA.Mass()))/(2.*M_PI)/ScalarConstant(std::string("vev")));
    
    kcpl0 = -M_I*ghgg;
    kcpl1 = kcpl0;
    
    vertex[vanz].cpl[0]  = kcpl0;
    vertex[vanz].cpl[1]  = kcpl0;
    vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
      
    vertex[vanz].Color.push_back(Color_Function(cf::G));     
    vertex[vanz].Color.back().SetParticleArg(0,2);     
    vertex[vanz].Color.back().SetStringArg('0','2');     
    vertex[vanz].Color.back().SetParticleArg(2,0);     
    vertex[vanz].Color.back().SetStringArg('2','0');     

    vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("PseudoTriangle",LF_Key()));     
    vertex[vanz].Lorentz.back()->SetParticleArg(0,2);     
      
    vertex[vanz].on      = 1;
    vertex.push_back(Single_Vertex());vanz++;
  }
}

void Interaction_Model_MSSM_EHC::c_SSS(std::vector<Single_Vertex>& vertex,int& vanz)  
{ 
  p_mothdm->c_SSS(vertex,vanz);
  p_moinos->c_SSS(vertex,vanz);
  p_moslepton->c_SSS(vertex,vanz);
  p_mosqcd->c_SSS(vertex,vanz);
  p_mosquark->c_SSS(vertex,vanz);
}

void Interaction_Model_MSSM_EHC::c_SSV(std::vector<Single_Vertex>& vertex,int& vanz) 
{ 
  p_mothdm->c_SSV(vertex,vanz);
  p_moinos->c_SSV(vertex,vanz);
  p_moslepton->c_SSV(vertex,vanz);
  p_mosqcd->c_SSV(vertex,vanz);
  p_mosquark->c_SSV(vertex,vanz);
}

void Interaction_Model_MSSM_EHC::c_SSVV(std::vector<Single_Vertex>& vertex,int& vanz) 
{ 
  p_mothdm->c_SSVV(vertex,vanz);
  p_moinos->c_SSVV(vertex,vanz);
  p_moslepton->c_SSVV(vertex,vanz);
  p_mosqcd->c_SSVV(vertex,vanz);
  p_mosquark->c_SSVV(vertex,vanz);
}

void Interaction_Model_MSSM_EHC::c_SSSS(std::vector<Single_Vertex>& vertex,int& vanz) 
{ 
  p_mothdm->c_SSSS(vertex,vanz);
  p_moinos->c_SSSS(vertex,vanz);
  p_moslepton->c_SSSS(vertex,vanz);
  p_mosqcd->c_SSSS(vertex,vanz);
  p_mosquark->c_SSSS(vertex,vanz);
  p_moslesqu->c_SSSS(vertex,vanz);
}

Interaction_Model_MSSM_EHC::~Interaction_Model_MSSM_EHC()
{
  delete p_mothdm;
  delete p_moinos;
  delete p_moslepton;
  delete p_mosqcd;
  delete p_mosquark;
  delete p_moslesqu;
}
