#include "Interaction_Model_SM_Phantom_U1.H"
#include "MathTools.H"
#include "Message.H"
#include "Run_Parameter.H"
#include <stdio.h>


using namespace AMEGIC;
using namespace ATOOLS;
using namespace std;

Interaction_Model_SM_Phantom_U1::Interaction_Model_SM_Phantom_U1(MODEL::Model_Base * _model,
					   string _cplscheme,string _yukscheme) :
  Interaction_Model_Base(_model,_cplscheme,_yukscheme)
{ 
  p_mosm  = new Interaction_Model_SM(p_model,_cplscheme,_yukscheme); 
  double Ecms2  = sqr(rpa.gen.Ecms());
  double hmass2 = sqr(Flavour(kf::h0).Mass());
  double Hmass2 = sqr(Flavour(kf::h0).Mass());

  g1    = Kabbala(string("g_1"),
		  sqrt(4.*M_PI*ScalarFunction(std::string("alpha_QED"),Ecms2)));
  g2    = Kabbala(string("g_1/\\sin\\theta_W"), 
		  g1.Value()/sqrt(ScalarConstant(std::string("sin2_thetaW"))));
  sintW = Kabbala(std::string("\\sin\\theta_W"),
		  sqrt(ScalarConstant(std::string("sin2_thetaW"))));
  costW = Kabbala(std::string("\\cos\\theta_W"),
		  sqrt(1.-ScalarConstant(std::string("sin2_thetaW"))));
  PL    = Kabbala(string("P_L"),1.);
  PR    = Kabbala(string("P_R"),1.);
  M_I   = Kabbala(string("i"),Complex(0.,1.));
  root2 = Kabbala(string("\\sqrt{2}"),sqrt(2.));
  vev   = Kabbala(string("v_{EW}"),ScalarConstant(string("vev")));
  GF    = Kabbala(string("GF"),ScalarConstant(string("GF")));
  tanb  = Kabbala(string("\\tan\\beta"),ScalarConstant(string("Tan(Beta)")));

  geffh = Kabbala(std::string("I_S^{(h)}"),
		  ScalarConstant(std::string("Higgs_GG_eff_h"))*
		  ScalarFunction(std::string("alpha_S"),hmass2)/
		  (2.*M_PI*ScalarConstant(std::string("vev"))));
  geffH = Kabbala(std::string("I_S^{(H)}"),
		  ScalarConstant(std::string("Higgs_GG_eff_H"))*
		  ScalarFunction(std::string("alpha_S"),Hmass2)/
		  (2.*M_PI*ScalarConstant(std::string("vev"))));
}

void Interaction_Model_SM_Phantom_U1::c_FFV(vector<Single_Vertex>& vertex,int& vanz)
{
  p_mosm->c_FFV(vertex,vanz);
}

void Interaction_Model_SM_Phantom_U1::c_VVV(vector<Single_Vertex>& vertex,int& vanz)
{
  p_mosm->c_VVV(vertex,vanz);
}
void Interaction_Model_SM_Phantom_U1::c_VVVV(vector<Single_Vertex>& vertex,int& vanz)
{
  p_mosm->c_VVVV(vertex,vanz);
}

void Interaction_Model_SM_Phantom_U1::c_FFS(vector<Single_Vertex>& vertex,int& vanz)  { 
  Flavour flh0(kf::h0), flH0(kf::H0);
  if (!flh0.IsOn() && !flH0.IsOn()) return;
  Kabbala kcpl0,kcpl1,massf,mixh,mixH;
  mixh = Kabbala(string("O_{11}"),ComplexMatrixElement("HiggsMix",0,0));
  mixH = Kabbala(string("O_{21}"),ComplexMatrixElement("HiggsMix",1,0));

  for (short int i=1;i<17;i++) {
    if (i==7) i=11;
    Flavour flav = Flavour(kf::code(i));
    if (flav.IsOn() && flav.IsFermion() && (flav.Yuk() > 0.)) {
      
      massf = Kabbala(string("M_{")+flav.TexName()+string("}(m_h^2)"),
		    ScalarFunction(string("m")+string(flav.Name()),sqr(flh0.Mass())));
      kcpl0 = -M_I*massf*mixh/vev;
      kcpl1 = kcpl0;
      if (!ATOOLS::IsZero(kcpl0.Value())) {
	vertex[vanz].in[0]   = flav;
	vertex[vanz].in[1]   = flh0;
	vertex[vanz].in[2]   = flav;
	vertex[vanz].cpl[0]  = kcpl0.Value();
	vertex[vanz].cpl[1]  = kcpl1.Value();
	vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;
	vertex[vanz].ncf     = 1;
	if (flav.Strong()) {
	  vertex[vanz].Color = new Color_Function(cf::D);     
	  vertex[vanz].Color->SetParticleArg(0,2);     
	  vertex[vanz].Color->SetStringArg('0','2');     
	}
	else 
	  vertex[vanz].Color  = new Color_Function(cf::None);	
	vertex[vanz].nlf      = 1;
	vertex[vanz].Lorentz  = new Lorentz_Function(lf::FFS);
	vertex[vanz].on       = 1;
	vertex.push_back(Single_Vertex());vanz++;
      }


      massf = Kabbala(string("M_{")+flav.TexName()+string("}(m_H^2)"),
		    ScalarFunction(string("m")+string(flav.Name()),sqr(flH0.Mass())));
      kcpl0 = -M_I*massf*mixH/vev;
      kcpl1 = kcpl0;
      if (!ATOOLS::IsZero(kcpl0.Value())) {
	vertex[vanz].in[0] = flav;
	vertex[vanz].in[1] = flH0;
	vertex[vanz].in[2] = flav;
	vertex[vanz].cpl[0]  = kcpl0.Value();
	vertex[vanz].cpl[1]  = kcpl1.Value();
	vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;
	vertex[vanz].ncf   = 1;
	if (flav.Strong()) {
	  vertex[vanz].Color = new Color_Function(cf::D);     
	  vertex[vanz].Color->SetParticleArg(0,2);     
	  vertex[vanz].Color->SetStringArg('0','2');     
	}
	else 
	  vertex[vanz].Color = new Color_Function(cf::None);
	vertex[vanz].nlf     = 1;
	vertex[vanz].Lorentz = new Lorentz_Function(lf::FFS);
	vertex[vanz].on      = 1;
	vertex.push_back(Single_Vertex());vanz++;
      }
    }
  }
}

void Interaction_Model_SM_Phantom_U1::c_VVS(vector<Single_Vertex>& vertex,int& vanz)  { 
  Flavour flh0(kf::h0), flH0(kf::H0);
  if (!flh0.IsOn() && !flH0.IsOn()) return;
  Kabbala kcpl0,kcpl1,massf,mixh,mixH;
  mixh = Kabbala(string("O_{11}"),ComplexMatrixElement("HiggsMix",0,0));
  mixH = Kabbala(string("O_{21}"),ComplexMatrixElement("HiggsMix",1,0));

  Kabbala num_2 = Kabbala(string("2"),2.);  
  Flavour flav(kf::W);
  // W h W
  if (flav.IsOn()) {
    vertex[vanz].in[0] = flav;
    vertex[vanz].in[1] = flh0;
    vertex[vanz].in[2] = flav;
    kcpl0 = M_I*g2*flav.Yuk()*mixh;
    kcpl1 = kcpl0;
    vertex[vanz].cpl[0]  = kcpl0.Value();
    vertex[vanz].cpl[1]  = vertex[vanz].cpl[0];
    vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
    vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;
    vertex[vanz].ncf   = 1;
    vertex[vanz].Color = new Color_Function(cf::None);     
    vertex[vanz].nlf     = 1;
    vertex[vanz].Lorentz = new Lorentz_Function(lf::Gab);     
    vertex[vanz].Lorentz->SetParticleArg(0,2);     
    vertex[vanz].on      = 1;
    vertex.push_back(Single_Vertex());vanz++;


    vertex[vanz].in[0] = flav;
    vertex[vanz].in[1] = flH0;
    vertex[vanz].in[2] = flav;
    kcpl0 = M_I*g2*flav.Yuk()*mixH;
    kcpl1 = kcpl0;
    vertex[vanz].cpl[0]  = kcpl0.Value();
    vertex[vanz].cpl[1]  = vertex[vanz].cpl[0];
    vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
    vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;
    vertex[vanz].ncf   = 1;
    vertex[vanz].Color = new Color_Function(cf::None);     
    vertex[vanz].nlf     = 1;
    vertex[vanz].Lorentz = new Lorentz_Function(lf::Gab);     
    vertex[vanz].Lorentz->SetParticleArg(0,2);     
    vertex[vanz].on      = 1;
    vertex.push_back(Single_Vertex());vanz++;
  }


  flav = Flavour(kf::Z);
  // Z h Z
  if (flav.IsOn()) {
    vertex[vanz].in[0] = flav;
    vertex[vanz].in[1] = flh0;
    vertex[vanz].in[2] = flav;
    kcpl0 = M_I*g2*flav.Yuk()/costW;
    kcpl1 = kcpl0;
    vertex[vanz].cpl[0]  = kcpl0.Value();
    vertex[vanz].cpl[1]  = kcpl1.Value();
    vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
    vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;
    vertex[vanz].ncf   = 1;
    vertex[vanz].Color = new Color_Function(cf::None);     
    vertex[vanz].nlf     = 1;
    vertex[vanz].Lorentz = new Lorentz_Function(lf::Gab);  
    vertex[vanz].Lorentz->SetParticleArg(0,2);     
    vertex[vanz].on      = 1;
    vertex.push_back(Single_Vertex());vanz++;


    vertex[vanz].in[0] = flav;
    vertex[vanz].in[1] = flH0;
    vertex[vanz].in[2] = flav;
    kcpl0 = M_I*g2*flav.Yuk()*mixH/costW;
    kcpl1 = kcpl0;
    vertex[vanz].cpl[0]  = kcpl0.Value();
    vertex[vanz].cpl[1]  = kcpl1.Value();
    vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
    vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;
    vertex[vanz].ncf   = 1;
    vertex[vanz].Color = new Color_Function(cf::None);     
    vertex[vanz].nlf     = 1;
    vertex[vanz].Lorentz = new Lorentz_Function(lf::Gab);  
    vertex[vanz].Lorentz->SetParticleArg(0,2);     
    vertex[vanz].on      = 1;
    vertex.push_back(Single_Vertex());vanz++;
  }


  flav = Flavour(kf::photon);
  // Photon h Photon
  if (flav.IsOn()) {
    vertex[vanz].in[0] = flav;
    vertex[vanz].in[1] = flh0;
    vertex[vanz].in[2] = flav;
    kcpl0 = M_I*geffh*mixh;
    kcpl1 = kcpl0;
    vertex[vanz].cpl[0]  = kcpl0.Value();
    vertex[vanz].cpl[1]  = kcpl0.Value();
    vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
    vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;
    vertex[vanz].ncf   = 1;
    vertex[vanz].Color = new Color_Function(cf::None);     
    vertex[vanz].nlf     = 1;
    vertex[vanz].Lorentz = new Lorentz_Function(lf::Triangle);     
    vertex[vanz].Lorentz->SetParticleArg(0,2);     
    vertex[vanz].on      = 1;
    vertex.push_back(Single_Vertex());vanz++;


    vertex[vanz].in[0] = flav;
    vertex[vanz].in[1] = flH0;
    vertex[vanz].in[2] = flav;    
    kcpl0 = M_I*geffH*mixH;
    kcpl1 = kcpl0;
    vertex[vanz].cpl[0]  = kcpl0.Value();
    vertex[vanz].cpl[1]  = kcpl0.Value();
    vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
    vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;
    vertex[vanz].ncf   = 1;
    vertex[vanz].Color = new Color_Function(cf::None);     
    vertex[vanz].nlf     = 1;
    vertex[vanz].Lorentz = new Lorentz_Function(lf::Triangle);     
    vertex[vanz].Lorentz->SetParticleArg(0,2);     
    vertex[vanz].on      = 1;
    vertex.push_back(Single_Vertex());vanz++;
  }

  Flavour flg(kf::gluon);
  // Gluon h Gluon
  if (flg.IsOn()) {
    vertex[vanz].in[0] = flg;
    vertex[vanz].in[1] = flh0;
    vertex[vanz].in[2] = flg;
    kcpl0 = M_I*geffh*mixh;
    kcpl1 = kcpl0;
    vertex[vanz].cpl[0]  = kcpl0.Value();
    vertex[vanz].cpl[1]  = kcpl0.Value();
    vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
    vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;
    vertex[vanz].ncf   = 1;
    vertex[vanz].Color = new Color_Function(cf::G);     
    vertex[vanz].Color->SetParticleArg(0,2);     
    vertex[vanz].Color->SetStringArg('0','2');     
    vertex[vanz].nlf     = 1;
    vertex[vanz].Lorentz = new Lorentz_Function(lf::Triangle);     
    vertex[vanz].Lorentz->SetParticleArg(0,2);     
    vertex[vanz].on      = 1;
    vertex.push_back(Single_Vertex());
    vanz++;

    vertex[vanz].in[0] = flg;
    vertex[vanz].in[1] = flH0;
    vertex[vanz].in[2] = flg;
    kcpl0 = M_I*geffH*mixH;
    kcpl1 = kcpl0;
    vertex[vanz].cpl[0]  = kcpl0.Value();
    vertex[vanz].cpl[1]  = kcpl0.Value();
    vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
    vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;
    vertex[vanz].ncf   = 1;
    vertex[vanz].Color = new Color_Function(cf::G);     
    vertex[vanz].Color->SetParticleArg(0,2);     
    vertex[vanz].Color->SetStringArg('0','2');     
    vertex[vanz].nlf     = 1;
    vertex[vanz].Lorentz = new Lorentz_Function(lf::Triangle);     
    vertex[vanz].Lorentz->SetParticleArg(0,2);     
    vertex[vanz].on      = 1;
    vertex.push_back(Single_Vertex());
    vanz++;
  }


  Flavour flsh(kf::shgluon);
  // gluon h shgluon
  if (flg.IsOn() && flsh.IsOn()) {
    vertex[vanz].in[2] = flg;
    vertex[vanz].in[0] = flsh;
    vertex[vanz].in[1] = flh0;    
    kcpl0 = M_I*mixh;
    kcpl1 = kcpl0;
    vertex[vanz].cpl[0]  = kcpl0.Value();
    vertex[vanz].cpl[1]  = kcpl0.Value();
    vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
    vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;
    vertex[vanz].ncf   = 1;
    vertex[vanz].Color = new Color_Function(cf::G);     
    vertex[vanz].Color->SetParticleArg(0,2);     
    vertex[vanz].Color->SetStringArg('0','2');     
    vertex[vanz].nlf     = 1;
    vertex[vanz].Lorentz = new Lorentz_Function(lf::C4GS);     
    vertex[vanz].Lorentz->SetParticleArg(0,2);     
    vertex[vanz].on      = 1;
    vertex[vanz].t       = -1;
    vertex.push_back(Single_Vertex());vanz++;


    vertex[vanz].in[2] = flg;
    vertex[vanz].in[0] = flsh;
    vertex[vanz].in[1] = flH0;    
    kcpl0 = M_I*mixH;
    kcpl1 = kcpl0;
    vertex[vanz].cpl[0]  = kcpl0.Value();
    vertex[vanz].cpl[1]  = kcpl0.Value();
    vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
    vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;
    vertex[vanz].ncf   = 1;
    vertex[vanz].Color = new Color_Function(cf::G);     
    vertex[vanz].Color->SetParticleArg(0,2);     
    vertex[vanz].Color->SetStringArg('0','2');     
    vertex[vanz].nlf     = 1;
    vertex[vanz].Lorentz = new Lorentz_Function(lf::C4GS);     
    vertex[vanz].Lorentz->SetParticleArg(0,2);     
    vertex[vanz].on      = 1;
    vertex[vanz].t       = -1;
    vertex.push_back(Single_Vertex());vanz++;
  }
}

void Interaction_Model_SM_Phantom_U1::c_SSS(vector<Single_Vertex>& vertex,int& vanz)  { 
  Flavour flh0(Flavour(kf::h0)), flH0(Flavour(kf::H0)), flA0(Flavour(kf::A0));

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
  massh2  = Kabbala(string("m_h^2"),sqr(flh0.Mass()));
  massH2  = Kabbala(string("m_H^2"),sqr(flH0.Mass()));
  num_2   = Kabbala(string("2"),2.);
  num_3   = Kabbala(string("3"),3.);

  vertex[vanz].in[0] = flh0;
  vertex[vanz].in[1] = flA0;
  vertex[vanz].in[2] = flA0;
  kcpl0              = -M_I*tanb*mix12*(massh2/vev);
  kcpl1              = kcpl0;
  vertex[vanz].cpl[0]  = kcpl0.Value();
  vertex[vanz].cpl[1]  = kcpl1.Value();
  vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
  vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;
  vertex[vanz].ncf     = 1;
  vertex[vanz].Color   = new Color_Function(cf::None);     
  vertex[vanz].nlf     = 1;
  vertex[vanz].Lorentz = new Lorentz_Function(lf::SSS);     
  vertex[vanz].on      = 1;
  vertex.push_back(Single_Vertex());vanz++;

  vertex[vanz].in[0] = flH0;
  vertex[vanz].in[1] = flA0;
  vertex[vanz].in[2] = flA0;
  kcpl0              = -M_I*tanb*mix22*(massH2/vev);
  kcpl1              = kcpl0;
  vertex[vanz].cpl[0]  = kcpl0.Value();
  vertex[vanz].cpl[1]  = kcpl1.Value();
  vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
  vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;
  vertex[vanz].ncf     = 1;
  vertex[vanz].Color   = new Color_Function(cf::None);     
  vertex[vanz].nlf     = 1;
  vertex[vanz].Lorentz = new Lorentz_Function(lf::SSS);     
  vertex[vanz].on      = 1;
  vertex.push_back(Single_Vertex());vanz++;
  
  vertex[vanz].in[0] = flh0;
  vertex[vanz].in[1] = flh0;
  vertex[vanz].in[2] = flH0;
  kcpl0              = M_I*(num_2*massh2+massH2)/vev*mix11*mix12*(mix11+mix21*tanb);
  kcpl1              = kcpl0;
  vertex[vanz].cpl[0]  = kcpl0.Value();
  vertex[vanz].cpl[1]  = kcpl1.Value();
  vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
  vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;
  vertex[vanz].ncf     = 1;
  vertex[vanz].Color   = new Color_Function(cf::None);     
  vertex[vanz].nlf     = 1;
  vertex[vanz].Lorentz = new Lorentz_Function(lf::SSS);     
  vertex[vanz].on      = 1;
  vertex.push_back(Single_Vertex());vanz++;
  // std::cout<<"-------------------------------"<<kcpl0.Value()
  // 	   <<" = "<<(abs(((num_2*massh2+massH2)/vev).Value()))<<"*"
  // 	   <<sqr(abs((mix11*mix12*(mix11+mix21*tanb)).Value()))<<endl
  // 	   <<"   from "<<massH2.Value()<<"/"<<massh2.Value()<<"/"<<(abs(vev.Value()))<<endl;


  vertex[vanz].in[0] = flH0;
  vertex[vanz].in[1] = flH0;
  vertex[vanz].in[2] = flh0;
  kcpl0              = M_I*(num_2*massH2+massh2)/vev*mix21*mix22*(mix12+mix22*tanb);
  kcpl1              = kcpl0;
  vertex[vanz].cpl[0]  = kcpl0.Value();
  vertex[vanz].cpl[1]  = kcpl1.Value();
  vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
  vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;
  vertex[vanz].ncf     = 1;
  vertex[vanz].Color   = new Color_Function(cf::None);     
  vertex[vanz].nlf     = 1;
  vertex[vanz].Lorentz = new Lorentz_Function(lf::SSS);     
  vertex[vanz].on      = 1;
  vertex.push_back(Single_Vertex());vanz++;


  vertex[vanz].in[0] = flh0;
  vertex[vanz].in[1] = flh0;
  vertex[vanz].in[2] = flh0;
  kcpl0              = -(num_3*M_I*massh2/vev)*(tanb*mix12_3+mix11_3);
  kcpl1              = kcpl0;
  vertex[vanz].cpl[0]  = kcpl0.Value();
  vertex[vanz].cpl[1]  = kcpl1.Value();
  vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
  vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;
  vertex[vanz].ncf     = 1;
  vertex[vanz].Color   = new Color_Function(cf::None);     
  vertex[vanz].nlf     = 1;
  vertex[vanz].Lorentz = new Lorentz_Function(lf::SSS);     
  vertex[vanz].on      = 1;
  vertex.push_back(Single_Vertex());vanz++;


  vertex[vanz].in[0] = flH0;
  vertex[vanz].in[1] = flH0;
  vertex[vanz].in[2] = flH0;
  kcpl0              = -(num_3*M_I*massH2/vev)*(tanb*mix22_3+mix21_3);
  kcpl1              = kcpl0;
  vertex[vanz].cpl[0]  = kcpl0.Value();
  vertex[vanz].cpl[1]  = kcpl1.Value();
  vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
  vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;
  vertex[vanz].ncf     = 1;
  vertex[vanz].Color   = new Color_Function(cf::None);     
  vertex[vanz].nlf     = 1;
  vertex[vanz].Lorentz = new Lorentz_Function(lf::SSS);     
  vertex[vanz].on      = 1;
  vertex.push_back(Single_Vertex());vanz++;
}

void Interaction_Model_SM_Phantom_U1::c_SSVV(vector<Single_Vertex>& vertex,int& vanz) { 
}

void Interaction_Model_SM_Phantom_U1::c_SSSS(vector<Single_Vertex>& vertex,int& vanz) { 
}

Interaction_Model_SM_Phantom_U1::~Interaction_Model_SM_Phantom_U1()
{
  delete p_mosm;
}
