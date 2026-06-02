#include "MODEL/DarkQCD/Model.H"
#include "ATOOLS/Phys/KF_Table.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "MODEL/Main/Running_AlphaQED.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "MODEL/Main/Strong_Coupling.H"
#include "MODEL/Main/Running_Fermion_Mass.H"
#include "MODEL/Main/Single_Vertex.H"
#include "PDF/Main/ISR_Handler.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Run_Parameter.H"

using namespace MODEL;
using namespace ATOOLS;
using namespace std;

DECLARE_GETTER(MODEL::DarkQCD,"DarkQCD",MODEL::Model_Base,MODEL::Model_Arguments);

Model_Base *ATOOLS::Getter<MODEL::Model_Base,MODEL::Model_Arguments,MODEL::DarkQCD>::
operator()(const Model_Arguments &args) const
{
  return new DarkQCD();
}



void ATOOLS::Getter<MODEL::Model_Base,MODEL::Model_Arguments,MODEL::DarkQCD>::
PrintInfo(ostream &str,const size_t width) const
{
  str<<"The Standard Model + dark SU(N)\n";
  str<<setw(width+4)<<" "<<"{\n"
     <<setw(width+7)<<" "<<"# possible parameters in yaml configuration [usage: \"keyword: value\"]\n"
     <<setw(width+7)<<" "<<"- EW_SCHEME (EW input scheme, see documentation)\n"
     <<setw(width+7)<<" "<<"- EW_REN_SCHEME (EW renormalisation scheme, see documentation)\n"
     <<setw(width+7)<<" "<<"- WIDTH_SCHEME (Fixed or CMS, see documentation)\n"
     <<setw(width+7)<<" "<<"- ALPHAS(MZ) (strong coupling at MZ)\n"
     <<setw(width+7)<<" "<<"- ORDER_ALPHAS (0,1,2 -> 1, 2, 3-loop running)\n"
     <<setw(width+7)<<" "<<"- 1/ALPHAQED(0) (alpha QED Thompson limit)\n"
     <<setw(width+7)<<" "<<"- ALPHAQED_DEFAULT_SCALE (scale for alpha_QED default)\n"
     <<setw(width+7)<<" "<<"- SIN2THETAW (weak mixing angle)\n"
     <<setw(width+7)<<" "<<"- VEV (Higgs vev)\n"
     <<setw(width+7)<<" "<<"- CKM_ORDER (0,1,2,3 - order of CKM expansion in Cabibbo angle)\n"
     <<setw(width+7)<<" "<<"- CKM_CABIBBO (Cabibbo angle in Wolfenstein parameterization)\n"
     <<setw(width+7)<<" "<<"- CKM_A (Wolfenstein A)\n"
     <<setw(width+7)<<" "<<"- CKM_RHO (Wolfenstein Rho)\n"
     <<setw(width+7)<<" "<<"- CKM_ETA (Wolfenstein Eta)\n"
     <<setw(width+7)<<" "<<"- CKM_ELEMENT[<i>][<j>] (explicit value for element, supersedes parametrisation)\n"
     <<setw(width+4)<<" "<<"}";
  str<<"Infrared continuation of alphaS:\n";
  str<<setw(width+4)<<" "<<"{\n"
     <<setw(width+7)<<" "<<"- AS_FORM (values 0,1,2,3,10, see documentation)\n"
     <<setw(width+7)<<" "<<"- Q2_AS (corresponding infrared parameter, see documentation)\n"
     <<setw(width+4)<<" "<<"}";
}

DarkQCD::DarkQCD() :
  Standard_Model() {
  m_name="DarkQCD";
  ParticleZprimeInit();
  AddParticle(kf_jv,0.,0.,0.,0,1, 2,1,1,1,0,"jv","jv","jv","jv",1,1);
  s_kftable[kf_jv]->Clear();
  s_kftable[kf_jv]->Add(Flavour(kf_gv));
  s_kftable[kf_jv]->Add(Flavour(kf_qv));
  s_kftable[kf_jv]->Add(Flavour(kf_qv).Bar());
  s_kftable[kf_jv]->Add(Flavour(kf_Dv));
  s_kftable[kf_jv]->Add(Flavour(kf_Dv).Bar());
  AddParticle(kf_jjv,0.,0.,0.,0,1, 2,1,1,1,0,"jjv","jjv","jjv","jjv",1,1);
  Settings& s = Settings::GetMainSettings();
  const double jet_mass_threshold{ s["JET_MASS_THRESHOLD"].Get<double>() };
  for (int i=1;i<7;i++) {
    Flavour addit((kf_code)i);
    if ((addit.Mass()==0.0 || !addit.IsMassive()) && addit.IsOn()) {
      if (addit.Mass(true)<=jet_mass_threshold) {
        s_kftable[kf_jjv]->Add(addit);
        s_kftable[kf_jjv]->Add(addit.Bar());
      }
      else {
        msg_Info()<<"Ignoring "<<addit<<" due to JJV_MASS_THRESHOLD.\n";
      }
    }
  }
  s_kftable[kf_jjv]->Add(Flavour(kf_gluon));
  s_kftable[kf_jjv]->Add(Flavour(kf_gv));
  s_kftable[kf_jjv]->Add(Flavour(kf_qv));
  s_kftable[kf_jjv]->Add(Flavour(kf_qv).Bar());
  s_kftable[kf_jjv]->Add(Flavour(kf_Dv));
  s_kftable[kf_jjv]->Add(Flavour(kf_Dv).Bar());
}

bool DarkQCD::ModelInit() {
  bool ret = Standard_Model::ModelInit();
  const Scoped_Settings& s = Settings::GetMainSettings()["Dark_QCD"];
  int    order_alphaS   = s["ORDER_ALPHAS"].SetDefault(2).Get<int>();
  int    th_alphaS      = s["THRESHOLD_ALPHAS"].SetDefault(1).Get<int>();
  double alphaDark = s["ALPHAS(MZp)"].SetDefault(0.118).Get<double>();
  double MZp = Flavour(kf_Zp).Mass();
  Running_AlphaDark * adark = new Running_AlphaDark(alphaDark,MZp,order_alphaS,th_alphaS,*p_isrhandlermap);
  p_constants->insert(make_pair(string("alpha_Dark"),alphaDark));
  p_functions->insert(make_pair(string("alpha_Dark"),adark));
  return ret;
}

void DarkQCD::InitVertices() {
  Standard_Model::InitVertices();
  InitZprimeVertices();
  InitDarkQCDVertices();
}


void DarkQCD::ParticleZprimeInit() {
  // add Zprime
  // kf_code,mass,radius,width,3*charge,strong,spin,majorana,take,stable,massive,idname,antiname,texname,antitexname
  AddParticle(kf_Zp,1000.,0,10.,0,0,2,-1,1,0,1,"Zprime","Zprime","Z^{\\prime}","Z^{\\prime}");
  AddParticle(kf_Dv,10,.0,.0,-1,3,1,0,1,1,0,"Dv","Dvb", "Dvd", "\\bar{Dv}");
  AddParticle(kf_qv,10,.0,.0,-1,3,1,0,1,1,0,"qv","qvb", "qvd", "\\bar{qv}");
  AddParticle(kf_gv,.0,.0,.0,0,8,2,-1,1,1,0,"gv","gv", "gv", "gv");
}


void DarkQCD::InitZprimeVertices() {
  if (!Flavour(kf_Zp).IsOn()) return;
  Kabbala g1("g_1_dark",sqrt(4.*M_PI*1/100));
  Kabbala cpl=g1*Kabbala("i",Complex(0.,1.));
  for (short int i=1;i<6;++i) {
    Flavour flav((kf_code)i);
    if (flav.IsOn() && flav.Charge()) {
      Kabbala Q("Q_dark_{"+flav.TexName()+"}",flav.Charge());
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle(flav.Bar());
      m_v.back().AddParticle(flav);
      m_v.back().AddParticle(Flavour(kf_Zp));
      m_v.back().Color.push_back(Color_Function(cf::D,1,2));
      m_v.back().Lorentz.push_back("FFV");
      m_v.back().cpl.push_back(cpl*Q);
      m_v.back().order[1]=1;
    }
  }

  Flavour flav(kf_qv);
  if(flav.IsOn()) {
    Kabbala Q("Qdark_{"+flav.TexName()+"}",flav.Charge());
    m_v.push_back(Single_Vertex());
    m_v.back().AddParticle(flav.Bar());
    m_v.back().AddParticle(flav);
    m_v.back().AddParticle(Flavour(kf_Zp));
    m_v.back().Color.push_back(Color_Function(cf::D,1,2));
    m_v.back().Lorentz.push_back("FFV");
    m_v.back().cpl.push_back(cpl*Q);
    m_v.back().order[1]=1;
  }
}

void DarkQCD::InitDarkQCDVertices() {
  Settings& s = Settings::GetMainSettings();
  if (!Flavour(kf_gv).IsOn()) return;
  m_dec_dark_g4 = false; //s["DECOMPOSE_4G_VERTEX"].Get<int>();
  Kabbala g3("g_3_dark",sqrt(4.*M_PI*ScalarConstant("alpha_Dark")));
  Kabbala cpl0=g3*Kabbala("i",Complex(0.,1.));
  Flavour flav;
  flav = Flavour(kf_qv);
  if (flav.IsOn()) {
    m_v.push_back(Single_Vertex());
    m_v.back().AddParticle(flav.Bar());
    m_v.back().AddParticle(flav);
    m_v.back().AddParticle(Flavour(kf_gv));
    m_v.back().Color.push_back(Color_Function(cf::T,3,2,1));
    m_v.back().Lorentz.push_back("FFV");
    m_v.back().cpl.push_back(cpl0);
    /* m_v.back().order[0]=1; */
    m_v.back().order.push_back(1);
  }
  flav = Flavour(kf_Dv);
  if (flav.IsOn()) {
    m_v.push_back(Single_Vertex());
    m_v.back().AddParticle(flav.Bar());
    m_v.back().AddParticle(flav);
    m_v.back().AddParticle(Flavour(kf_gv));
    m_v.back().Color.push_back(Color_Function(cf::T,3,2,1));
    m_v.back().Lorentz.push_back("FFV");
    m_v.back().cpl.push_back(cpl0);
    /* m_v.back().order[0]=1; */
    m_v.back().order.push_back(1);
  }
  Kabbala cpl1=-g3;
  m_v.push_back(Single_Vertex());
  for (size_t i(0);i<3;++i) m_v.back().AddParticle(Flavour(kf_gv));
  m_v.back().Color.push_back(Color_Function(cf::F,1,2,3));
  m_v.back().Lorentz.push_back("VVV");
  m_v.back().cpl.push_back(cpl1);
  m_v.back().order[0]=1;
  // if (m_dec_dark_g4) {
  //   m_v.push_back(Single_Vertex());
  //   for (size_t i(0);i<2;++i) m_v.back().AddParticle(Flavour(kf_gv));
  //   m_v.back().AddParticle(Flavour(kf_gv_qgc));
  //   m_v.back().Color.push_back(Color_Function(cf::F,1,2,3));
  //   m_v.back().Lorentz.push_back("VVP");
  //   m_v.back().cpl.push_back(cpl1);
  //   m_v.back().order[0]=1;
  //   m_v.back().dec=1;
  // }
  Kabbala cpl2=g3*g3*Kabbala("i",Complex(0.,1.));
  m_v.push_back(Single_Vertex());
  for (size_t i(0);i<4;++i) m_v.back().AddParticle(Flavour(kf_gv));
  for (size_t i(0);i<3;++i) m_v.back().cpl.push_back(cpl2);
  m_v.back().Color.push_back
    (Color_Function(cf::F,-1,1,2,new Color_Function(cf::F,3,4,-1)));
  m_v.back().Color.push_back
    (Color_Function(cf::F,-1,1,3,new Color_Function(cf::F,2,4,-1)));
  m_v.back().Color.push_back
    (Color_Function(cf::F,-1,1,4,new Color_Function(cf::F,2,3,-1)));
  m_v.back().Lorentz.push_back("VVVVA");
  m_v.back().Lorentz.push_back("VVVVB");
  m_v.back().Lorentz.push_back("VVVVC");
  /* m_v.back().order[0]=2; */
  m_v.back().order.push_back(1);
  if (m_dec_dark_g4) m_v.back().dec=-1;

}
