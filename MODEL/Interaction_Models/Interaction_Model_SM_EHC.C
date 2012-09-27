#include "MODEL/Interaction_Models/Interaction_Model_SM_EHC.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Data_Reader.H"

using namespace MODEL;
using namespace ATOOLS;
using namespace std;

DECLARE_GETTER(Interaction_Model_SM_EHC_Getter,"SM+EHC",
	       Interaction_Model_Base,Interaction_Model_Arguments);

Interaction_Model_Base *Interaction_Model_SM_EHC_Getter::
operator()(const Interaction_Model_Arguments &args) const
{
  return new Interaction_Model_SM_EHC
    (args.p_model,args.m_cplscheme,args.m_yukscheme);
}

void Interaction_Model_SM_EHC_Getter::PrintInfo
(std::ostream &str,const size_t width) const
{ 
  str<<"The Standard Model + Effective Higgs Coupling"; 
}

Interaction_Model_SM_EHC::Interaction_Model_SM_EHC(MODEL::Model_Base * _model,
					       std::string _cplscheme,std::string _yukscheme) :
  Interaction_Model_Base("SM+EHC",_model,_cplscheme,_yukscheme)
{ 
  m_loops=true;
  p_mosm      = new Interaction_Model_SM(p_model,_cplscheme,_yukscheme); 
  
  double scale = rpa->gen.CplScale();
  
  g1    = Kabbala(string("g_1"),
		  sqrt(4.*M_PI*ScalarFunction(std::string("alpha_QED"),scale)));
  if(ScalarNumber(std::string("WidthScheme"))==0){
    g2       = Kabbala(string("g_1/\\sin\\theta_W"), 
		     g1.Value()/sqrt(ScalarConstant(std::string("sin2_thetaW"))));
    sintW    = Kabbala(std::string("\\sin\\theta_W"),
		     sqrt(ScalarConstant(std::string("sin2_thetaW"))));
    costW    = Kabbala(std::string("\\cos\\theta_W"),
		     sqrt(1.-ScalarConstant(std::string("sin2_thetaW"))));
    vev      = Kabbala(string("v_{EW}"),ScalarConstant(std::string("vev")));
  }else{
    g2       = Kabbala(string("g_1/\\sin\\theta_W"), 
		     g1.Value()/sqrt(ComplexConstant(std::string("csin2_thetaW"))));
    sintW    = Kabbala(std::string("\\sin\\theta_W"),
		     sqrt(ComplexConstant(std::string("csin2_thetaW"))));
    costW    = Kabbala(std::string("\\cos\\theta_W"),
		     sqrt(1.-ComplexConstant(std::string("csin2_thetaW"))));
    vev      = Kabbala(string("v_{EW}"),ComplexConstant(std::string("cvev")));
  }
  
  Data_Reader read(" ",";","!","=");
  double ehc_scale2 = read.GetValue<double>("EHC_SCALE2", 
					    sqr(Flavour(kf_h0).Mass()));
  p_model->GetScalarConstants()->insert(make_pair(string("EHC_SCALE2"),
						  ehc_scale2));
  // h photon photon coupling
  double aqedpph=read.GetValue<double>("1/ALPHAQED_PPH",
				       ScalarFunction(string("alpha_QED"),
						      ehc_scale2));
  ghpp  = Kabbala(std::string("ghpp"),ScalarConstant(string("h0_pp_fac"))*
                  aqedpph/(2.*M_PI*vev.Value()));
  // h g g coupling
  //msg_Out()<<"Evaluate asggh now with |"
  //	   <<ScalarFunction(string("alpha_S"),ehc_scale2)<<"|.\n";
  double asggh=read.GetValue<double>("ALPHAS_GGH",
				     ScalarFunction(string("alpha_S"),
						    ehc_scale2));
  //msg_Out()<<" ... yields asggH = "<<asggh<<".\n";
  ghgg  = Kabbala(std::string("ghgg"),ScalarConstant(string("h0_gg_fac"))*
		  asggh/(2.*M_PI*vev.Value()));
  //msg_Info()<<METHOD<<"() {\n"
  //    <<"  ggh coupling is "<<ghgg.Value()<<" ( \\alpha_s = "<<asggh<<" )\n"
  //    <<"  pph coupling is "<<ghpp.Value()<<" ( 1/\\alpha_qed = "
  //<<1./aqedpph<<" )\n"
  //    <<"  taken at \\mu = "<<sqrt(ehc_scale2)<<"\n}\n";
  g3  = Kabbala(string("g_3"),
		sqrt(4.*M_PI*ScalarFunction(std::string("alpha_S"),scale)));
  PL    = Kabbala(string("P_L"),1.);
  PR    = Kabbala(string("P_R"),1.);
  M_I   = Kabbala(string("i"),Complex(0.,1.));
  root2 = Kabbala(string("\\sqrt{2}"),sqrt(2.));
}

void Interaction_Model_SM_EHC::c_FFV(std::vector<Single_Vertex>& vertex,int& vanz)
{
  p_mosm->c_FFV(vertex,vanz);
}

void Interaction_Model_SM_EHC::c_VVV(std::vector<Single_Vertex>& vertex,int& vanz)
{
  p_mosm->c_VVV(vertex,vanz);
}
void Interaction_Model_SM_EHC::c_VVVV(std::vector<Single_Vertex>& vertex,int& vanz)
{
  p_mosm->c_VVVV(vertex,vanz);
  
  Kabbala kcpl0,kcpl1; 
  Flavour flh(kf_h0);
  Flavour flg(kf_gluon);
  if(!flh.IsOn()||!flg.IsOn())return;
  
  // 3 gluon higgs
  vertex[vanz].nleg    = 4;
  for (short int i=0;i<3;i++) vertex[vanz].in[i] = flg;
  vertex[vanz].in[3] = flh;
  
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

  Flavour flsh(kf_shgluon);
  kcpl0 = M_I*g3*g3*ghgg; 
  kcpl1 = kcpl0; 
  if(!flsh.IsOn()) return;
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

void Interaction_Model_SM_EHC::c_FFS(std::vector<Single_Vertex>& vertex,int& vanz) 
{
  p_mosm->c_FFS(vertex,vanz);
}
void Interaction_Model_SM_EHC::c_VVS(std::vector<Single_Vertex>& vertex,int& vanz) 
{
  p_mosm->c_VVS(vertex,vanz);
  
  Flavour flh(kf_h0);
  Kabbala kcpl0,kcpl1;  
  Kabbala num_2 = Kabbala(string("2"),2.);  
  
  if (!flh.IsOn()) return;
  
  Flavour flav(kf_photon);
  // Photon h Photon
  if (flav.IsOn()) {
    vertex[vanz].in[0] = flav;
    vertex[vanz].in[1] = flh;
    vertex[vanz].in[2] = flav;
    
    kcpl0 = M_I*ghpp;
    kcpl1 = kcpl0;
    
    vertex[vanz].cpl[0]  = kcpl0;
    vertex[vanz].cpl[1]  = kcpl0;
    vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
    
    vertex[vanz].Color.push_back(Color_Function(cf::None));     
    
    vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("Triangle",LF_Key()));     
    vertex[vanz].Lorentz.back()->SetParticleArg(0,2);     
    
    vertex[vanz].on      = 1;
    vertex.push_back(Single_Vertex());vanz++;
  }
  
  Flavour flg(kf_gluon);
  // Gluon h Gluon
  if (flg.IsOn()) {
    vertex[vanz].in[0] = flg;
    vertex[vanz].in[1] = flh;
    vertex[vanz].in[2] = flg;
    
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
void Interaction_Model_SM_EHC::c_SSS(std::vector<Single_Vertex>& vertex,int& vanz) 
{
  p_mosm->c_SSS(vertex,vanz);
}
void Interaction_Model_SM_EHC::c_SSSS(std::vector<Single_Vertex>& vertex,int& vanz) 
{ 
  p_mosm->c_SSSS(vertex,vanz); 
}
void Interaction_Model_SM_EHC::c_SSVV(std::vector<Single_Vertex>& vertex,int& vanz) 
{ 
  p_mosm->c_SSVV(vertex,vanz); 
}


Interaction_Model_SM_EHC::~Interaction_Model_SM_EHC()
{
  delete  p_mosm;
}

















