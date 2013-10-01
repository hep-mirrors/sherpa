#include "MODEL/Interaction_Models/Interaction_Model_SM_Phantom_U1.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include <stdio.h>


using namespace MODEL;
using namespace ATOOLS;
using namespace std;

DECLARE_GETTER(Interaction_Model_SM_Phantom_U1,"SM+Phantom_U1",
	       Interaction_Model_Base,Interaction_Model_Arguments);

Interaction_Model_Base *ATOOLS::Getter
<Interaction_Model_Base,Interaction_Model_Arguments,Interaction_Model_SM_Phantom_U1>::
operator()(const Interaction_Model_Arguments &args) const
{
  return new Interaction_Model_SM_Phantom_U1
    (args.p_model,args.m_cplscheme,args.m_yukscheme);
}

void ATOOLS::Getter<Interaction_Model_Base,Interaction_Model_Arguments,
		    Interaction_Model_SM_Phantom_U1>::PrintInfo
(std::ostream &str,const size_t width) const
{ 
  str<<"The Standard Model + phantom Higgs"; 
}

Interaction_Model_SM_Phantom_U1::
Interaction_Model_SM_Phantom_U1(MODEL::Model_Base * _model,
				string _cplscheme,string _yukscheme) :
  Interaction_Model_Base("SM+Phantom_U1",_model,_cplscheme,_yukscheme)
{ 
  m_loops=true;
  p_moew  = new Interaction_Model_EW(p_model,_cplscheme,_yukscheme); 
  p_moqcd = new Interaction_Model_QCD(p_model,_cplscheme,_yukscheme); 
  double hmass2 = sqr(Flavour(kf_h0).Mass());
  double Hmass2 = sqr(Flavour(kf_H0).Mass());

  g1    = Kabbala(string("g_1"),
		  sqrt(4.*M_PI*ScalarFunction(std::string("alpha_QED"),rpa->gen.CplScale())));
  if(ScalarNumber(std::string("WidthScheme"))==0){
    g2       = Kabbala(string("g_1/\\sin\\theta_W"), 
  		     g1.Value()/sqrt(ScalarConstant(std::string("sin2_thetaW"))));
    sintW    = Kabbala(std::string("\\sin\\theta_W"),
  		     sqrt(ScalarConstant(std::string("sin2_thetaW"))));
    costW    = Kabbala(std::string("\\cos\\theta_W"),
  		     sqrt(1.-ScalarConstant(std::string("sin2_thetaW"))));
    vev      = Kabbala(string("v_{Higgs_SM}"),ScalarConstant(std::string("vev")));
  }else{
    g2       = Kabbala(string("g_1/\\sin\\theta_W"), 
		     g1.Value()/sqrt(ComplexConstant(std::string("csin2_thetaW"))));
    sintW    = Kabbala(std::string("\\sin\\theta_W"),
		     sqrt(ComplexConstant(std::string("csin2_thetaW"))));
    costW    = Kabbala(std::string("\\cos\\theta_W"),
		     sqrt(1.-ComplexConstant(std::string("csin2_thetaW"))));
    vev      = Kabbala(string("v_{Higgs_SM}"),ComplexConstant(std::string("cvev")));
  }
  PL    = Kabbala(string("P_L"),1.);
  PR    = Kabbala(string("P_R"),1.);
  M_I   = Kabbala(string("i"),Complex(0.,1.));
  root2 = Kabbala(string("\\sqrt{2}"),sqrt(2.));
  tanb  = Kabbala(string("\\tan\\beta"),ScalarConstant(string("Tan(Beta)")));

  ghgg = Kabbala(std::string("I_S^{(h)}"),
		 ComplexConstant(std::string("h0_gg_fac"))*
		 ScalarFunction(std::string("alpha_S"),hmass2)/
		 (2.*M_PI*vev.Value()));
  gHgg = Kabbala(std::string("I_S^{(H)}"),
		 ComplexConstant(std::string("H0_gg_fac"))*
		 ScalarFunction(std::string("alpha_S"),Hmass2)/
		 (2.*M_PI*vev.Value()));
  ghpp = Kabbala(std::string("I_P^{(h)}"),
		 ComplexConstant(std::string("h0_pp_fac"))*
		 ScalarFunction(std::string("alpha_QED"),hmass2)/
		 (2.*M_PI*vev.Value()));
  gHpp = Kabbala(std::string("I_P^{(H)}"),
		 ComplexConstant(std::string("H0_pp_fac"))*
		 ScalarFunction(std::string("alpha_QED"),Hmass2)/
		 (2.*M_PI*vev.Value()));
}

void Interaction_Model_SM_Phantom_U1::c_FFV(vector<Single_Vertex>& vertex,int& vanz)
{
  p_moew->c_FFV(vertex,vanz);
  p_moqcd->c_FFV(vertex,vanz);
}

void Interaction_Model_SM_Phantom_U1::c_VVV(vector<Single_Vertex>& vertex,int& vanz)
{
  p_moew->c_VVV(vertex,vanz);
  p_moqcd->c_VVV(vertex,vanz);
}
void Interaction_Model_SM_Phantom_U1::c_VVVV(vector<Single_Vertex>& vertex,int& vanz)
{
  p_moew->c_VVVV(vertex,vanz);
  p_moqcd->c_VVVV(vertex,vanz);
}

void Interaction_Model_SM_Phantom_U1::c_FFS(vector<Single_Vertex>& vertex,int& vanz)  { 
  Flavour flh0(kf_h0), flH0(kf_H0);
  if (!flh0.IsOn() && !flH0.IsOn()) return;
  Kabbala kcpl0,kcpl1,massf,mixh,mixH;
  mixh = Kabbala(string("O_{11}"),ComplexMatrixElement("HiggsMix",0,0));
  mixH = Kabbala(string("O_{21}"),ComplexMatrixElement("HiggsMix",1,0));

  for (short int i=1;i<17;i++) {
    if (i==7) i=11;
    Flavour flav = Flavour((kf_code)(i));
    if (flav.IsOn() && flav.IsFermion() && (flav.Yuk() > 0.)) {
      massf = Kabbala(string("M_{")+flav.TexName()+string("}(m_h^2)"),
		      sqrt(sqr(ScalarFunction(string("m")+flav.IDName(),sqr(flh0.Mass())))-
		      Complex(0.,1.)*flav.Width()*ScalarFunction(string("m")+flav.IDName(),sqr(flh0.Mass()))));

      kcpl0 = -M_I*massf*mixh/vev;
      kcpl1 = kcpl0;
      if (!ATOOLS::IsZero(kcpl0.Value())) {
	vertex[vanz].in[0]   = flav;
	vertex[vanz].in[1]   = flh0;
	vertex[vanz].in[2]   = flav;
	vertex[vanz].cpl[0]  = kcpl0;
	vertex[vanz].cpl[1]  = kcpl1;
	vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	
	if (flav.Strong()) {
	  vertex[vanz].Color.push_back(Color_Function(cf::D));     
	  vertex[vanz].Color.back().SetParticleArg(0,2);     
	  vertex[vanz].Color.back().SetStringArg('0','2');     
	}
	else 
	  vertex[vanz].Color.push_back(Color_Function(cf::None));	
	
	vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("FFS",LF_Key()));
	vertex[vanz].on       = 1;
	vertex.push_back(Single_Vertex());vanz++;
      }

      kcpl0 = -M_I*massf*mixH/vev;
      kcpl1 = kcpl0;
      if (!ATOOLS::IsZero(kcpl0.Value())) {
	vertex[vanz].in[0] = flav;
	vertex[vanz].in[1] = flH0;
	vertex[vanz].in[2] = flav;
	vertex[vanz].cpl[0]  = kcpl0;
	vertex[vanz].cpl[1]  = kcpl1;
	vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	
	if (flav.Strong()) {
	  vertex[vanz].Color.push_back(Color_Function(cf::D));     
	  vertex[vanz].Color.back().SetParticleArg(0,2);     
	  vertex[vanz].Color.back().SetStringArg('0','2');     
	}
	else 
	  vertex[vanz].Color.push_back(Color_Function(cf::None));
	
	vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("FFS",LF_Key()));
	vertex[vanz].on      = 1;
	vertex.push_back(Single_Vertex());vanz++;
      }
    }
  }
}

void Interaction_Model_SM_Phantom_U1::c_VVS(vector<Single_Vertex>& vertex,int& vanz)  { 
  Flavour flh0(kf_h0), flH0(kf_H0);
  if (!flh0.IsOn() && !flH0.IsOn()) return;
  Kabbala kcpl0,kcpl1,massf,mixh,mixH;
  mixh = Kabbala(string("O_{11}"),ComplexMatrixElement("HiggsMix",0,0));
  mixH = Kabbala(string("O_{21}"),ComplexMatrixElement("HiggsMix",1,0));

  Kabbala num_2 = Kabbala(string("2"),2.);  
  Flavour flav(kf_Wplus);
  // W h W
  if (flav.IsOn()) {
    vertex[vanz].in[0] = flav.Bar();
    vertex[vanz].in[1] = flh0;
    vertex[vanz].in[2] = flav.Bar();
    if(ScalarNumber(std::string("WidthScheme"))==0){
      kcpl0 = M_I*g2*flav.Yuk()*mixh;
    }else{
      kcpl0 = M_I*g2*sqrt(sqr(flav.Yuk())-Complex(0.,1.)*flav.Width()*flav.Yuk())*mixh;
    }
    kcpl1 = kcpl0;
    vertex[vanz].cpl[0]  = kcpl0;
    vertex[vanz].cpl[1]  = vertex[vanz].cpl[0];
    vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
    
    vertex[vanz].Color.push_back(Color_Function(cf::None));     
    
    vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("Gab",LF_Key()));     
    vertex[vanz].Lorentz.back()->SetParticleArg(0,2);     
    vertex[vanz].on      = 1;
    vertex.push_back(Single_Vertex());vanz++;


    vertex[vanz].in[0] = flav;
    vertex[vanz].in[1] = flH0;
    vertex[vanz].in[2] = flav;
    if(ScalarNumber(std::string("WidthScheme"))==0){
      kcpl0 = M_I*g2*flav.Yuk()*mixH;
    }else{
      kcpl0 = M_I*g2*sqrt(sqr(flav.Yuk())-Complex(0.,1.)*flav.Width()*flav.Yuk())*mixH;
    }
    kcpl1 = kcpl0;
    vertex[vanz].cpl[0]  = kcpl0;
    vertex[vanz].cpl[1]  = vertex[vanz].cpl[0];
    vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
    
    vertex[vanz].Color.push_back(Color_Function(cf::None));     
    
    vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("Gab",LF_Key()));     
    vertex[vanz].Lorentz.back()->SetParticleArg(0,2);     
    vertex[vanz].on      = 1;
    vertex.push_back(Single_Vertex());vanz++;
  }


  flav = Flavour(kf_Z);
  // Z h Z
  if (flav.IsOn()) {
    vertex[vanz].in[0] = flav;
    vertex[vanz].in[1] = flh0;
    vertex[vanz].in[2] = flav;
    if(ScalarNumber(std::string("WidthScheme"))==0){
      kcpl0 = M_I*g2*flav.Yuk()/costW;
    }else{
      kcpl0 = M_I*g2*sqrt(sqr(flav.Yuk())-Complex(0.,1.)*flav.Width()*flav.Yuk())/costW;
    }
    kcpl1 = kcpl0;
    vertex[vanz].cpl[0]  = kcpl0;
    vertex[vanz].cpl[1]  = kcpl1;
    vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
    
    vertex[vanz].Color.push_back(Color_Function(cf::None));     
    
    vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("Gab",LF_Key()));  
    vertex[vanz].Lorentz.back()->SetParticleArg(0,2);     
    vertex[vanz].on      = 1;
    vertex.push_back(Single_Vertex());vanz++;


    vertex[vanz].in[0] = flav;
    vertex[vanz].in[1] = flH0;
    vertex[vanz].in[2] = flav;
    if(ScalarNumber(std::string("WidthScheme"))==0){
      kcpl0 = M_I*g2*flav.Yuk()*mixH/costW;
    }else{
      kcpl0 = M_I*g2*sqrt(sqr(flav.Yuk())-Complex(0.,1.)*flav.Width()*flav.Yuk())*mixH/costW;
    }
    kcpl1 = kcpl0;
    vertex[vanz].cpl[0]  = kcpl0;
    vertex[vanz].cpl[1]  = kcpl1;
    vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
    
    vertex[vanz].Color.push_back(Color_Function(cf::None));     
    
    vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("Gab",LF_Key()));  
    vertex[vanz].Lorentz.back()->SetParticleArg(0,2);     
    vertex[vanz].on      = 1;
    vertex.push_back(Single_Vertex());vanz++;
  }

  bool ehc(true);

  flav = Flavour(kf_photon);
  // Photon h Photon
  if (flav.IsOn()) {
    vertex[vanz].in[0] = flav;
    vertex[vanz].in[1] = flh0;
    vertex[vanz].in[2] = flav;
    kcpl0 = M_I*ghpp*mixh;
    kcpl1 = kcpl0;
    vertex[vanz].cpl[0]  = kcpl0;
    vertex[vanz].cpl[1]  = kcpl0;
    vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
    
    vertex[vanz].Color.push_back(Color_Function(cf::None));     
    
    vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("Triangle",LF_Key()));     
    vertex[vanz].Lorentz.back()->SetParticleArg(0,2);     
    vertex[vanz].on      = ehc;
    vertex.push_back(Single_Vertex());vanz++;


    vertex[vanz].in[0] = flav;
    vertex[vanz].in[1] = flH0;
    vertex[vanz].in[2] = flav;    
    kcpl0 = M_I*gHpp*mixH;
    kcpl1 = kcpl0;
    vertex[vanz].cpl[0]  = kcpl0;
    vertex[vanz].cpl[1]  = kcpl0;
    vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
    
    vertex[vanz].Color.push_back(Color_Function(cf::None));     
    
    vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("Triangle",LF_Key()));     
    vertex[vanz].Lorentz.back()->SetParticleArg(0,2);     
    vertex[vanz].on      = ehc;
    vertex.push_back(Single_Vertex());vanz++;
  }

  ehc = true;
  Flavour flg(kf_gluon);
  // Gluon h Gluon
  if (flg.IsOn()) {
    vertex[vanz].in[0] = flg;
    vertex[vanz].in[1] = flh0;
    vertex[vanz].in[2] = flg;
    kcpl0 = M_I*ghgg*mixh;
    kcpl1 = kcpl0;
    vertex[vanz].cpl[0]  = kcpl0;
    vertex[vanz].cpl[1]  = kcpl0;
    vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
    
    vertex[vanz].Color.push_back(Color_Function(cf::G));     
    vertex[vanz].Color.back().SetParticleArg(0,2);     
    vertex[vanz].Color.back().SetStringArg('0','2');     
    
    vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("Triangle",LF_Key()));     
    vertex[vanz].Lorentz.back()->SetParticleArg(0,2);     
    vertex[vanz].on      = ehc;
    vertex.push_back(Single_Vertex());
    vanz++;

    vertex[vanz].in[0] = flg;
    vertex[vanz].in[1] = flH0;
    vertex[vanz].in[2] = flg;
    kcpl0 = M_I*gHgg*mixH;
    kcpl1 = kcpl0;
    vertex[vanz].cpl[0]  = kcpl0;
    vertex[vanz].cpl[1]  = kcpl0;
    vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
    
    vertex[vanz].Color.push_back(Color_Function(cf::G));     
    vertex[vanz].Color.back().SetParticleArg(0,2);     
    vertex[vanz].Color.back().SetStringArg('0','2');     
    
    vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("Triangle",LF_Key()));     
    vertex[vanz].Lorentz.back()->SetParticleArg(0,2);     
    vertex[vanz].on      = ehc;
    vertex.push_back(Single_Vertex());
    vanz++;
  }


  Flavour flsh(kf_shgluon);
  // gluon h shgluon
  if (flg.IsOn() && flsh.IsOn()) {
    vertex[vanz].in[2] = flg;
    vertex[vanz].in[0] = flsh;
    vertex[vanz].in[1] = flh0;    
    kcpl0 = M_I*mixh;
    kcpl1 = kcpl0;
    vertex[vanz].cpl[0]  = kcpl0;
    vertex[vanz].cpl[1]  = kcpl0;
    vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
    
    vertex[vanz].Color.push_back(Color_Function(cf::G));     
    vertex[vanz].Color.back().SetParticleArg(0,2);     
    vertex[vanz].Color.back().SetStringArg('0','2');     
    
    vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("C4GS",LF_Key()));     
    vertex[vanz].Lorentz.back()->SetParticleArg(0,2);     
    vertex[vanz].on      = ehc;
    vertex[vanz].t       = -1;
    vertex.push_back(Single_Vertex());vanz++;


    vertex[vanz].in[2] = flg;
    vertex[vanz].in[0] = flsh;
    vertex[vanz].in[1] = flH0;    
    kcpl0 = M_I*mixH;
    kcpl1 = kcpl0;
    vertex[vanz].cpl[0]  = kcpl0;
    vertex[vanz].cpl[1]  = kcpl0;
    vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
    
    vertex[vanz].Color.push_back(Color_Function(cf::G));     
    vertex[vanz].Color.back().SetParticleArg(0,2);     
    vertex[vanz].Color.back().SetStringArg('0','2');     
    
    vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("C4GS",LF_Key()));     
    vertex[vanz].Lorentz.back()->SetParticleArg(0,2);     
    vertex[vanz].on      = ehc;
    vertex[vanz].t       = -1;
    vertex.push_back(Single_Vertex());vanz++;
  }
}

void Interaction_Model_SM_Phantom_U1::c_SSS(vector<Single_Vertex>& vertex,int& vanz)  { 
  Flavour flh0(kf_h0), flH0(kf_H0), flA0(kf_A0);

  Kabbala kcpl0,kcpl1,massh2,massH2,mix11,mix21,mix12,mix22;
  Kabbala mix11_3,mix21_3,mix12_3,mix22_3,num_2,num_3;
  mix11   = Kabbala(string("O_{11}"),ComplexMatrixElement("HiggsMix",0,0));
  mix21   = Kabbala(string("O_{21}"),ComplexMatrixElement("HiggsMix",1,0));
  mix12   = Kabbala(string("O_{12}"),ComplexMatrixElement("HiggsMix",0,1));
  mix22   = Kabbala(string("O_{22}"),ComplexMatrixElement("HiggsMix",1,1));
  mix11_3 = Kabbala(string("O_{11}^3"),pow(ComplexMatrixElement("HiggsMix",0,0),3));
  mix21_3 = Kabbala(string("O_{21}^3"),pow(ComplexMatrixElement("HiggsMix",1,0),3));
  mix12_3 = Kabbala(string("O_{12}^3"),pow(ComplexMatrixElement("HiggsMix",0,1),3));
  mix22_3 = Kabbala(string("O_{22}^3"),pow(ComplexMatrixElement("HiggsMix",1,1),3));
  if(ScalarNumber(std::string("WidthScheme"))==0){
    massh2  = Kabbala(string("m_h^2"),sqr(flh0.Mass()));
    massH2  = Kabbala(string("m_H^2"),sqr(flH0.Mass()));
  }else{
    massh2  = Kabbala(string("m_h^2"),sqr(flh0.Mass())-Complex(0.,1.)*flh0.Mass()*flh0.Width());
    massH2  = Kabbala(string("m_H^2"),sqr(flH0.Mass())-Complex(0.,1.)*flH0.Mass()*flH0.Width());
  }
  num_2   = Kabbala(string("2"),2.);
  num_3   = Kabbala(string("3"),3.);

  vertex[vanz].in[0] = flh0;
  vertex[vanz].in[1] = flA0;
  vertex[vanz].in[2] = flA0;
  kcpl0              = -M_I*tanb*mix12*(massh2/vev);
  kcpl1              = kcpl0;
  vertex[vanz].cpl[0]  = kcpl0;
  vertex[vanz].cpl[1]  = kcpl1;
  vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
  
  vertex[vanz].Color.push_back(Color_Function(cf::None));     
  
  vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("SSS",LF_Key()));     
  vertex[vanz].on      = 1;
  vertex.push_back(Single_Vertex());vanz++;

  vertex[vanz].in[0] = flH0;
  vertex[vanz].in[1] = flA0;
  vertex[vanz].in[2] = flA0;
  kcpl0              = -M_I*tanb*mix22*(massH2/vev);
  kcpl1              = kcpl0;
  vertex[vanz].cpl[0]  = kcpl0;
  vertex[vanz].cpl[1]  = kcpl1;
  vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
  
  vertex[vanz].Color.push_back(Color_Function(cf::None));     
  
  vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("SSS",LF_Key()));     
  vertex[vanz].on      = 1;
  vertex.push_back(Single_Vertex());vanz++;
  
  vertex[vanz].in[0] = flh0;
  vertex[vanz].in[1] = flh0;
  vertex[vanz].in[2] = flH0;
  kcpl0              = M_I*(num_2*massh2+massH2)/vev*mix11*mix12*(mix11+mix21*tanb);
  kcpl1              = kcpl0;
  vertex[vanz].cpl[0]  = kcpl0;
  vertex[vanz].cpl[1]  = kcpl1;
  vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
  
  vertex[vanz].Color.push_back(Color_Function(cf::None));     
  
  vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("SSS",LF_Key()));     
  vertex[vanz].on      = 1;
  vertex.push_back(Single_Vertex());vanz++;


  vertex[vanz].in[0] = flH0;
  vertex[vanz].in[1] = flH0;
  vertex[vanz].in[2] = flh0;
  kcpl0              = M_I*(num_2*massH2+massh2)/vev*mix21*mix22*(mix12+mix22*tanb);
  kcpl1              = kcpl0;
  vertex[vanz].cpl[0]  = kcpl0;
  vertex[vanz].cpl[1]  = kcpl1;
  vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
  
  vertex[vanz].Color.push_back(Color_Function(cf::None));     
  
  vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("SSS",LF_Key()));     
  vertex[vanz].on      = 1;
  vertex.push_back(Single_Vertex());vanz++;


  vertex[vanz].in[0] = flh0;
  vertex[vanz].in[1] = flh0;
  vertex[vanz].in[2] = flh0;
  kcpl0              = -(num_3*M_I*massh2/vev)*(tanb*mix12_3+mix11_3);
  kcpl1              = kcpl0;
  vertex[vanz].cpl[0]  = kcpl0;
  vertex[vanz].cpl[1]  = kcpl1;
  vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
  
  vertex[vanz].Color.push_back(Color_Function(cf::None));     
  
  vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("SSS",LF_Key()));     
  vertex[vanz].on      = 1;
  vertex.push_back(Single_Vertex());vanz++;


  vertex[vanz].in[0] = flH0;
  vertex[vanz].in[1] = flH0;
  vertex[vanz].in[2] = flH0;
  kcpl0              = -(num_3*M_I*massH2/vev)*(tanb*mix22_3+mix21_3);
  kcpl1              = kcpl0;
  vertex[vanz].cpl[0]  = kcpl0;
  vertex[vanz].cpl[1]  = kcpl1;
  vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
  
  vertex[vanz].Color.push_back(Color_Function(cf::None));     
  
  vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("SSS",LF_Key()));     
  vertex[vanz].on      = 1;
  vertex.push_back(Single_Vertex());vanz++;
}

void Interaction_Model_SM_Phantom_U1::c_SSVV(vector<Single_Vertex>& vertex,int& vanz) { 
  Flavour flh0(kf_h0), flH0(kf_H0);
  Flavour flavW(kf_Wplus);
  Flavour flavZ(kf_Z);

  Kabbala kcpl0,kcpl1,massh2,massH2,mix11,mix21,mix12,mix22;
  Kabbala num_2,num_4;
  mix11   = Kabbala(string("O_{11}"),ComplexMatrixElement("HiggsMix",0,0));
  mix21   = Kabbala(string("O_{21}"),ComplexMatrixElement("HiggsMix",1,0));
  mix12   = Kabbala(string("O_{12}"),ComplexMatrixElement("HiggsMix",0,1));
  mix22   = Kabbala(string("O_{22}"),ComplexMatrixElement("HiggsMix",1,1));
  if(ScalarNumber(std::string("WidthScheme"))==0){
    massh2  = Kabbala(string("m_h^2"),sqr(flh0.Mass()));
    massH2  = Kabbala(string("m_H^2"),sqr(flH0.Mass()));
  }else{
    massh2  = Kabbala(string("m_h^2"),sqr(flh0.Mass())-Complex(0.,1.)*flh0.Mass()*flh0.Width());
    massH2  = Kabbala(string("m_H^2"),sqr(flH0.Mass())-Complex(0.,1.)*flH0.Mass()*flH0.Width());
  }
  num_2   = Kabbala(string("2"),2.); 
  num_4   = Kabbala(string("4"),4.);
 
  // h0 - Z - Z - h0   
  if (flavZ.IsOn() && flh0.IsOn()) {
    vertex[vanz].in[0] = flavZ;
    vertex[vanz].in[1] = flh0;
    vertex[vanz].in[2] = flh0;
    vertex[vanz].in[3] = flavZ;
    
    vertex[vanz].nleg     = 4;
    
    kcpl0 = (M_I*g2*g2/(costW*costW*num_2))*mix11*mix11;
    kcpl1 = kcpl0;
    
    vertex[vanz].cpl[0]  = kcpl0;
    vertex[vanz].cpl[1]  = kcpl1;
    vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
    
    
    vertex[vanz].Color.push_back(Color_Function(cf::None));     
    
    
    vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("VVSS",LF_Key()));
    vertex[vanz].Lorentz.back()->SetParticleArg(0,3);     
    
    vertex[vanz].on      = 1;
    vertex[vanz].oew     = 2;
    vertex.push_back(Single_Vertex());vanz++;
  }

 // H0 - Z - Z - h0 
  if (flavZ.IsOn() && flh0.IsOn() && flH0.IsOn()) {
    vertex[vanz].in[0] = flavZ;
    vertex[vanz].in[1] = flh0;
    vertex[vanz].in[2] = flH0;
    vertex[vanz].in[3] = flavZ;
    
    vertex[vanz].nleg     = 4;
    
    kcpl0 = (M_I*g2*g2/(costW*costW*num_2))*mix21*mix11;
    kcpl1 = kcpl0;
    
    vertex[vanz].cpl[0]  = kcpl0;
    vertex[vanz].cpl[1]  = kcpl1;
    vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
    
    
    vertex[vanz].Color.push_back(Color_Function(cf::None));     
    
    
    vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("VVSS",LF_Key()));
    vertex[vanz].Lorentz.back()->SetParticleArg(0,3);     
    
    vertex[vanz].on      = 1;
    vertex[vanz].oew     = 2;
    vertex.push_back(Single_Vertex());vanz++;
  }

 // H0 - Z - Z - H0   
  if (flavZ.IsOn() && flH0.IsOn()) {
    vertex[vanz].in[0] = flavZ;
    vertex[vanz].in[1] = flH0;
    vertex[vanz].in[2] = flH0;
    vertex[vanz].in[3] = flavZ;
    
    vertex[vanz].nleg     = 4;
    
    kcpl0 = (M_I*g2*g2/(costW*costW*num_2))*mix21*mix21;
    kcpl1 = kcpl0;
    
    vertex[vanz].cpl[0]  = kcpl0;
    vertex[vanz].cpl[1]  = kcpl1;
    vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
    
    
    vertex[vanz].Color.push_back(Color_Function(cf::None));     
    
    
    vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("VVSS",LF_Key()));
    vertex[vanz].Lorentz.back()->SetParticleArg(0,3);     
    
    vertex[vanz].on      = 1;
    vertex[vanz].oew     = 2;
    vertex.push_back(Single_Vertex());vanz++;
  }
  // h0 - W - W - h0  
  if (flavW.IsOn() && flh0.IsOn()) {
    vertex[vanz].in[0] = flavW.Bar();
    vertex[vanz].in[1] = flh0;
    vertex[vanz].in[2] = flh0;
    vertex[vanz].in[3] = flavW.Bar();
    
    vertex[vanz].nleg     = 4;
    
    kcpl0 = (M_I*g2*g2/num_2)*mix11*mix11;
    kcpl1 = kcpl0;
    
    vertex[vanz].cpl[0]  = kcpl0;
    vertex[vanz].cpl[1]  = kcpl1;
    vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
    
    
    vertex[vanz].Color.push_back(Color_Function(cf::None));     
    
    
    vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("VVSS",LF_Key()));     
    vertex[vanz].Lorentz.back()->SetParticleArg(0,3);     
    
    vertex[vanz].on      = 1;
    vertex[vanz].oew     = 2;
    vertex.push_back(Single_Vertex());vanz++;
  }
 // H0 - W - W - h0  
  if (flavW.IsOn() && flh0.IsOn() && flH0.IsOn()) {
    vertex[vanz].in[0] = flavW.Bar();
    vertex[vanz].in[1] = flH0;
    vertex[vanz].in[2] = flh0;
    vertex[vanz].in[3] = flavW.Bar();
    
    vertex[vanz].nleg     = 4;
    
    kcpl0 = (M_I*g2*g2/num_2)*mix21*mix11;
    kcpl1 = kcpl0;
    
    vertex[vanz].cpl[0]  = kcpl0;
    vertex[vanz].cpl[1]  = kcpl1;
    vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
    
    
    vertex[vanz].Color.push_back(Color_Function(cf::None));     
    
    
    vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("VVSS",LF_Key()));     
    vertex[vanz].Lorentz.back()->SetParticleArg(0,3);     
    
    vertex[vanz].on      = 1;
    vertex[vanz].oew     = 2;
    vertex.push_back(Single_Vertex());vanz++;
  }
 // H0 - W - W - H0  
  if (flavW.IsOn() && flH0.IsOn()) {
    vertex[vanz].in[0] = flavW.Bar();
    vertex[vanz].in[1] = flH0;
    vertex[vanz].in[2] = flH0;
    vertex[vanz].in[3] = flavW.Bar();
    
    vertex[vanz].nleg     = 4;
    
    kcpl0 = (M_I*g2*g2/num_2)*mix21*mix21;
    kcpl1 = kcpl0;
    
    vertex[vanz].cpl[0]  = kcpl0;
    vertex[vanz].cpl[1]  = kcpl1;
    vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
    
    
    vertex[vanz].Color.push_back(Color_Function(cf::None));     
    
    
    vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("VVSS",LF_Key()));     
    vertex[vanz].Lorentz.back()->SetParticleArg(0,3);     
    
    vertex[vanz].on      = 1;
    vertex[vanz].oew     = 2;
    vertex.push_back(Single_Vertex());vanz++;
  }
}

void Interaction_Model_SM_Phantom_U1::c_SSSS(vector<Single_Vertex>& vertex,int& vanz) { 
  Flavour flh0(kf_h0), flH0(kf_H0), flA0(kf_A0);
  Kabbala kcpl0,kcpl1,massh2,massH2,mix11,mix21,mix12,mix22;
  Kabbala num_2,num_4,num_3,num_6,num_8,num_24;
  Kabbala mix11_3,mix21_3,mix12_3,mix22_3;
  Kabbala mix11_4,mix21_4,mix12_4,mix22_4; 
  Kabbala mix11_2,mix21_2,mix12_2,mix22_2;
  mix11   = Kabbala(string("O_{11}"),ComplexMatrixElement("HiggsMix",0,0));
  mix21   = Kabbala(string("O_{21}"),ComplexMatrixElement("HiggsMix",1,0));
  mix12   = Kabbala(string("O_{12}"),ComplexMatrixElement("HiggsMix",0,1));
  mix22   = Kabbala(string("O_{22}"),ComplexMatrixElement("HiggsMix",1,1));
  mix11_3 = Kabbala(string("O_{11}^3"),pow(ComplexMatrixElement("HiggsMix",0,0),3));
  mix21_3 = Kabbala(string("O_{21}^3"),pow(ComplexMatrixElement("HiggsMix",1,0),3));
  mix12_3 = Kabbala(string("O_{12}^3"),pow(ComplexMatrixElement("HiggsMix",0,1),3));
  mix22_3 = Kabbala(string("O_{22}^3"),pow(ComplexMatrixElement("HiggsMix",1,1),3));				    
  mix11_4 = Kabbala(string("O_{11}^4"),pow(ComplexMatrixElement("HiggsMix",0,0),4));
  mix21_4 = Kabbala(string("O_{21}^4"),pow(ComplexMatrixElement("HiggsMix",1,0),4));
  mix12_4 = Kabbala(string("O_{12}^4"),pow(ComplexMatrixElement("HiggsMix",0,1),4));
  mix22_4 = Kabbala(string("O_{22}^4"),pow(ComplexMatrixElement("HiggsMix",1,1),4));	
  mix11_2 = Kabbala(string("O_{11}^2"),pow(ComplexMatrixElement("HiggsMix",0,0),2));
  mix21_2 = Kabbala(string("O_{21}^2"),pow(ComplexMatrixElement("HiggsMix",1,0),2));
  mix12_2 = Kabbala(string("O_{12}^2"),pow(ComplexMatrixElement("HiggsMix",0,1),2));
  mix22_2 = Kabbala(string("O_{22}^2"),pow(ComplexMatrixElement("HiggsMix",1,1),2));
  if(ScalarNumber(std::string("WidthScheme"))==0){
    massh2  = Kabbala(string("m_h^2"),sqr(flh0.Mass()));
    massH2  = Kabbala(string("m_H^2"),sqr(flH0.Mass()));
  }else{
    massh2  = Kabbala(string("m_h^2"),sqr(flh0.Mass())-Complex(0.,1.)*flh0.Mass()*flh0.Width());
    massH2  = Kabbala(string("m_H^2"),sqr(flH0.Mass())-Complex(0.,1.)*flH0.Mass()*flH0.Width());
  }
  num_2   = Kabbala(string("2"),2.); 
  num_3   = Kabbala(string("3"),3.);
  num_4   = Kabbala(string("4"),4.);
  num_8   = Kabbala(string("8"),8.);							    
  num_24  = Kabbala(string("24"),24.);
  num_6  = Kabbala(string("6"),6.);
  // h0-h0-h0-h0
  if (flh0.IsOn()) {  
    vertex[vanz].in[0] = flh0;
    vertex[vanz].in[1] = flh0;
    vertex[vanz].in[2] = flh0;
    vertex[vanz].in[3] = flh0;

    vertex[vanz].nleg  = 4;  
    
    kcpl0 = -M_I*num_3/(vev * vev)*
      (tanb*tanb*mix12_4*(mix11_2*massH2+mix12_2*massh2)
               + mix11_4*(mix11_2*massh2 + mix12_2*massH2) 
               - num_2*mix11_3*mix12_3*tanb*(massH2-massh2));
    kcpl1 = kcpl0;
    
    vertex[vanz].cpl[0]  = kcpl0;
    vertex[vanz].cpl[1]  = kcpl1;
    vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();

    
    vertex[vanz].Color.push_back(Color_Function(cf::None));     
    
    
    vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("SSSS",LF_Key()));     

    vertex[vanz].on      = 1;
    vertex[vanz].oew     = 2;
    vertex.push_back(Single_Vertex());vanz++;
  }
  // H0-H0-H0-H0
 if (flH0.IsOn()) {  
    vertex[vanz].in[0] = flH0;
    vertex[vanz].in[1] = flH0;
    vertex[vanz].in[2] = flH0;
    vertex[vanz].in[3] = flH0;

    vertex[vanz].nleg  = 4;  
    
    kcpl0 = -M_I*num_3/(vev * vev)*
      (tanb*tanb*mix22_4*(mix11_2*massH2+mix12_2*massh2)
               + mix21_4*(mix11_2*massh2 + mix12_2*massH2) 
               - num_2*mix11_3*mix12_3*tanb*(massH2-massh2));
    kcpl1 = kcpl0;
    
    vertex[vanz].cpl[0]  = kcpl0;
    vertex[vanz].cpl[1]  = kcpl1;
    vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();

    
    vertex[vanz].Color.push_back(Color_Function(cf::None));     
    
    
    vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("SSSS",LF_Key()));     

    vertex[vanz].on      = 1;
    vertex[vanz].oew     = 2;
    vertex.push_back(Single_Vertex());vanz++;
  }
  // H0-H0-H0-h0
 if (flH0.IsOn() && flh0.IsOn()) {  
    vertex[vanz].in[0] = flH0;
    vertex[vanz].in[1] = flH0;
    vertex[vanz].in[2] = flH0;
    vertex[vanz].in[3] = flh0;

    vertex[vanz].nleg  = 4;  
    
    kcpl0 = M_I*num_3*mix11*mix21/(vev * vev)*(mix12+mix11*tanb)*
            ((massH2*(mix22_3*tanb+mix21_3)+
              massh2*(mix21*mix11_2+mix22*mix12_2*tanb)));
    kcpl1 = kcpl0;
    
    vertex[vanz].cpl[0]  = kcpl0;
    vertex[vanz].cpl[1]  = kcpl1;
    vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();

    
    vertex[vanz].Color.push_back(Color_Function(cf::None));     
    
    
    vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("SSSS",LF_Key()));     

    vertex[vanz].on      = 1;
    vertex[vanz].oew     = 2;
    vertex.push_back(Single_Vertex());vanz++;
  }
 // h0-h0-h0-H0
 if (flH0.IsOn() && flh0.IsOn()) {  
    vertex[vanz].in[0] = flh0;
    vertex[vanz].in[1] = flh0;
    vertex[vanz].in[2] = flh0;
    vertex[vanz].in[3] = flH0;

    vertex[vanz].nleg  = 4;  
    
    kcpl0 = M_I*num_3*mix11*mix12/(vev * vev)*(mix22+mix21*tanb)*
      (massh2*(mix12_3*tanb+mix11_3)+
       massH2*(mix11*mix21_2+mix12*mix22_2*tanb));
    kcpl1 = kcpl0;
    
    vertex[vanz].cpl[0]  = kcpl0;
    vertex[vanz].cpl[1]  = kcpl1;
    vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();

    
    vertex[vanz].Color.push_back(Color_Function(cf::None));     
    
    
    vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("SSSS",LF_Key()));     

    vertex[vanz].on      = 1;
    vertex[vanz].oew     = 2;
    vertex.push_back(Single_Vertex());vanz++;
  }
 // h0-h0-H0-H0
 if (flH0.IsOn() && flh0.IsOn()) {  
    vertex[vanz].in[0] = flh0;
    vertex[vanz].in[1] = flh0;
    vertex[vanz].in[2] = flH0;
    vertex[vanz].in[3] = flH0;

    vertex[vanz].nleg  = 4;  
    
    kcpl0 = M_I*mix12*mix11/(vev * vev)*
            ( (massH2-massh2)*tanb*(mix11_4-num_4*mix11_2*mix12_2+mix12_4)
            - num_3*mix11*mix12*(massh2*mix11_2+massH2*mix12_2+
            (massh2*mix12_2+massH2*mix11_2)*tanb*tanb));
    kcpl1 = kcpl0;
    
    vertex[vanz].cpl[0]  = kcpl0;
    vertex[vanz].cpl[1]  = kcpl1;
    vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();

    
    vertex[vanz].Color.push_back(Color_Function(cf::None));     
    
    
    vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("SSSS",LF_Key()));     

    vertex[vanz].on      = 1;
    vertex[vanz].oew     = 2;
    vertex.push_back(Single_Vertex());vanz++;
  }
// h0-h0-A0-A0
 if (flA0.IsOn() && flh0.IsOn()) {  
    vertex[vanz].in[0] = flh0;
    vertex[vanz].in[1] = flh0;
    vertex[vanz].in[2] = flA0;
    vertex[vanz].in[3] = flA0;

    vertex[vanz].nleg  = 4;  
    
    kcpl0 = M_I/(vev * vev)*
            ( (massH2-massh2)*tanb*mix11*mix12*mix11*mix11-
              tanb*tanb*(massh2*mix12_2+massH2*mix11_2)*mix12*mix12);
    kcpl1 = kcpl0;
    
    vertex[vanz].cpl[0]  = kcpl0;
    vertex[vanz].cpl[1]  = kcpl1;
    vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();

    
    vertex[vanz].Color.push_back(Color_Function(cf::None));     
    
    
    vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("SSSS",LF_Key()));     

    vertex[vanz].on      = 1;
    vertex[vanz].oew     = 2;
    vertex.push_back(Single_Vertex());vanz++;
  }
// H0-H0-A0-A0
 if (flA0.IsOn() && flH0.IsOn()) {  
    vertex[vanz].in[0] = flH0;
    vertex[vanz].in[1] = flH0;
    vertex[vanz].in[2] = flA0;
    vertex[vanz].in[3] = flA0;

    vertex[vanz].nleg  = 4;  
    
    kcpl0 = M_I/(vev * vev)*
            ((massH2-massh2)*tanb*mix11*mix12*mix21*mix21-
              tanb*tanb*(massh2*mix12_2+massH2*mix11_2)*mix22*mix22);
    kcpl1 = kcpl0;
    
    vertex[vanz].cpl[0]  = kcpl0;
    vertex[vanz].cpl[1]  = kcpl1;
    vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();

    
    vertex[vanz].Color.push_back(Color_Function(cf::None));     
    
    
    vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("SSSS",LF_Key()));     

    vertex[vanz].on      = 1;
    vertex[vanz].oew     = 2;
    vertex.push_back(Single_Vertex());vanz++;
  }
// h0-H0-A0-A0
 if (flA0.IsOn() && flh0.IsOn() && flH0.IsOn()) {  
    vertex[vanz].in[0] = flh0;
    vertex[vanz].in[1] = flH0;
    vertex[vanz].in[2] = flA0;
    vertex[vanz].in[3] = flA0;

    vertex[vanz].nleg  = 4;  
    
    kcpl0 = M_I/(vev * vev)*
            ( (massH2-massh2)*tanb*mix11*mix12*mix11*mix21-
              tanb*tanb*(massh2*mix12_2+massH2*mix11_2)*mix12*mix22);
    kcpl1 = kcpl0;
    
    vertex[vanz].cpl[0]  = kcpl0;
    vertex[vanz].cpl[1]  = kcpl1;
    vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();

    
    vertex[vanz].Color.push_back(Color_Function(cf::None));     
    
    
    vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("SSSS",LF_Key()));     

    vertex[vanz].on      = 1;
    vertex[vanz].oew     = 2;
    vertex.push_back(Single_Vertex());vanz++;
  }
}

Interaction_Model_SM_Phantom_U1::~Interaction_Model_SM_Phantom_U1()
{
  delete p_moew;
  delete p_moqcd;
}
