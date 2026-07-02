#include "MODEL/DarkQCD/Model.H"
#include "ATOOLS/Phys/KF_Table.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "HADRONS++/Main/Hadron_Decay_Map.H"
#include "HADRONS++/Main/Hadron_Decay_Table.H"
#include "HADRONS++/Main/Hadron_Decay_Channel.H"
#include "HADRONS++/Main/Mixing_Handler.H"
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
  // AddParticle(kf_dark_j,10.,0.,0.,0,1, 2,1,1,1,0,"dark_j","dark_j","dark_j","dark_j",1,1);
  // s_kftable[kf_dark_j]->Clear();
  // s_kftable[kf_dark_j]->Add(Flavour(kf_dark_g));
  // s_kftable[kf_dark_j]->Add(Flavour(kf_dark_q1));
  // s_kftable[kf_dark_j]->Add(Flavour(kf_dark_q1).Bar());
  // s_kftable[kf_dark_j]->Add(Flavour(kf_dark_q2));
  // s_kftable[kf_dark_j]->Add(Flavour(kf_dark_q2).Bar());
  // //s_kftable[kf_dark_j]->Add(Flavour(kf_Dv));
  // //s_kftable[kf_dark_j]->Add(Flavour(kf_Dv).Bar());
  // AddParticle(kf_dark_jj,10.,0.,0.,0,1, 2,1,1,1,0,"dark_jj","dark_jj","dark_jj","dark_jj",1,1);
  // Settings& s = Settings::GetMainSettings();
  // const double jet_mass_threshold{ s["JET_MASS_THRESHOLD"].Get<double>() };
  // for (int i=1;i<7;i++) {
  //   Flavour addit((kf_code)i);
  //   if ((addit.Mass()==0.0 || !addit.IsMassive()) && addit.IsOn()) {
  //     if (addit.Mass(true)<=jet_mass_threshold) {
  //       s_kftable[kf_dark_jj]->Add(addit);
  //       s_kftable[kf_dark_jj]->Add(addit.Bar());
  //     }
  //     else {
  //       msg_Info()<<"Ignoring "<<addit<<" due to DARK_JJ_MASS_THRESHOLD.\n";
  //     }
  //   }
  // }
  // s_kftable[kf_dark_jj]->Add(Flavour(kf_gluon));
  // s_kftable[kf_dark_j]->Add(Flavour(kf_dark_g));
  // s_kftable[kf_dark_j]->Add(Flavour(kf_dark_q1));
  // s_kftable[kf_dark_j]->Add(Flavour(kf_dark_q1).Bar());
  // s_kftable[kf_dark_j]->Add(Flavour(kf_dark_q2));
  // s_kftable[kf_dark_j]->Add(Flavour(kf_dark_q2).Bar());
  //s_kftable[kf_dark_jj]->Add(Flavour(kf_Dv));
  //s_kftable[kf_dark_jj]->Add(Flavour(kf_Dv).Bar());
}

bool DarkQCD::ModelInit() {
  bool ret = Standard_Model::ModelInit();
  const Scoped_Settings& s = Settings::GetMainSettings()["Dark_QCD"];
  int    order_alphaS = s["ORDER_ALPHAS"].SetDefault(2).Get<int>();
  int    th_alphaS    = s["THRESHOLD_ALPHAS"].SetDefault(1).Get<int>();
  double alphaDark    = s["ALPHAS(MZp)"].SetDefault(0.118).Get<double>();
  double MZp          = Flavour(kf_Z0_2).Mass();
  Running_AlphaDark * adark = new Running_AlphaDark(alphaDark,MZp,order_alphaS,th_alphaS,*p_isrhandlermap);
  p_constants->insert(make_pair(string("alpha_Dark"),alphaDark));
  p_functions->insert(make_pair(string("alpha_Dark"),adark));
  m_fPi               = 0.137;
  m_gRhoPiPi          = 5.;
  m_fEta              = m_fPi/(16.*M_PI);
  m_gOmegaPiPi        = m_gRhoPiPi/(16.*M_PI);
  return ret;
}

void DarkQCD::InitVertices() {
  Standard_Model::InitVertices();
  InitZprimeVertices();
  InitDarkQCDVertices();
}


void DarkQCD::ParticleZprimeInit() {
  // add Zprime
  //          kf_code,       mass,radius,width,3*charge,strong,2*spin,
  //                     majorana,  take,stable,massive,idname,antiname;
  AddParticle(kf_Z0_2,      1000.,  0.0, 10.0,       0,    0,      2,
	                       -1, true, false,   true, "Zprime","Zprime","Z^{\\prime}","Z^{\\prime}");
  AddParticle(kf_dark_q1,      10,  0.0,  0.0,       0,    3,      1,
	                        0, true, true,    true, "dark_q1","dark_q1b", "\\tilde{q}_1", "\\bar{\\tilde{q}}_1");
  AddParticle(kf_dark_q2,      10,  0.0,  0.0,       0,    3,      1,
	                        0, true, true,    true, "dark_q2","dark_q2b", "\\tilde{q}_2", "\\bar{\\tilde{q}}_2");
  AddParticle(kf_dark_g,      0.0,  0.0,  0.0,       0,    8,      2,
	                       -1, true, true,   false,  "dark_g",  "dark_g",   "\\tilde{g}", "\\tilde{g}");
  //          kf_code,         mass,radius,width,3*charge,2*spin, take,stable,       idname,       antiname;
  AddParticle(kf_dark_pi11,     1.0,  0.0,  0.1,       0,    0,   true,     0,  "dark_pi11",   "dark_pi11");
  AddParticle(kf_dark_pi21,     1.1,  0.0,  0.0,       0,    0,   true,     1,  "dark_pi21",   "dark_pi21");
  AddParticle(kf_dark_eta22,    4.0,  0.0,  0.5,       0,    0,   true,     0,  "dark_eta22",  "dark_eta22");
  AddParticle(kf_dark_rho11,   18.0,  0.0,  5.0,       0,    2,   true,     0,  "dark_rho11",  "dark_rho11");
  AddParticle(kf_dark_rho21,   18.1,  0.0,  5.0,       0,    2,   true,     0,  "dark_rho21",  "dark_rho21");
  AddParticle(kf_dark_omega22, 22.0,  0.0,  3.0,       0,    2,   true,     0, "dark_omega22", "dark_omega22");
  if (Flavour(kf_dark_pi11).IsOn())    m_hadrons.insert(Flavour(kf_dark_pi11));
  if (Flavour(kf_dark_pi21).IsOn())    m_hadrons.insert(Flavour(kf_dark_pi21));
  if (Flavour(kf_dark_eta22).IsOn())   m_hadrons.insert(Flavour(kf_dark_eta22));
  if (Flavour(kf_dark_rho11).IsOn())   m_hadrons.insert(Flavour(kf_dark_rho11));
  if (Flavour(kf_dark_rho21).IsOn())   m_hadrons.insert(Flavour(kf_dark_rho21));
  if (Flavour(kf_dark_omega22).IsOn()) m_hadrons.insert(Flavour(kf_dark_omega22));
}


void DarkQCD::InitZprimeVertices() {
  if (!Flavour(kf_Z0_2).IsOn()) return;
  Kabbala g1("g_1_dark",sqrt(4.*M_PI*1./100.));
  Kabbala cpl=g1*Kabbala("i",Complex(0.,1.));
  for (short int i=1;i<=6;++i) {
    Flavour flav((kf_code)i);
    if (flav.IsOn() && flav.Charge() && Flavour(kf_Z0_2).IsOn()) {
      msg_Out()<<"Init Z'qq vertex for "<<flav<<"\n";
      Kabbala Q("Q_dark_{"+flav.TexName()+"}",flav.Charge());
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle(flav.Bar());
      m_v.back().AddParticle(flav);
      m_v.back().AddParticle(Flavour(kf_Z0_2));
      m_v.back().Color.push_back(Color_Function(cf::D,1,2));
      m_v.back().Lorentz.push_back("FFV");
      m_v.back().cpl.push_back(cpl);
      m_v.back().order[1]=1;
      msg_Out()<<m_v.back()<<"\n";
    }
  }
  for (short int i=11;i<=16;++i) {
    Flavour flav((kf_code)i);
    if (flav.IsOn() && flav.Charge() && Flavour(kf_Z0_2).IsOn()) {
      msg_Out()<<"Init Z'qq vertex for "<<flav<<"\n";
      Kabbala Q("Q_dark_{"+flav.TexName()+"}",flav.Charge());
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle(flav.Bar());
      m_v.back().AddParticle(flav);
      m_v.back().AddParticle(Flavour(kf_Z0_2));
      m_v.back().Color.push_back(Color_Function(cf::D,1,2));
      m_v.back().Lorentz.push_back("FFV");
      m_v.back().cpl.push_back(cpl);
      m_v.back().order[1]=1;
      msg_Out()<<m_v.back()<<"\n";
    }
  }

  for (long int i=4900101;i<4900103;++i) {
    Flavour flav((kf_code)i);
    if (flav.IsOn() && Flavour(kf_Z0_2).IsOn()) {
      msg_Out()<<"Init Z'qq vertex for "<<flav<<"\n";
      Kabbala Q("Q_dark_{"+flav.TexName()+"}",1.);
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle(flav.Bar());
      m_v.back().AddParticle(flav);
      m_v.back().AddParticle(Flavour(kf_Z0_2));
      m_v.back().Color.push_back(Color_Function(cf::D,1,2));
      m_v.back().Lorentz.push_back("FFV");
      m_v.back().cpl.push_back(cpl);
      m_v.back().order[1]=1;
      msg_Out()<<m_v.back()<<"\n";
    }
  }
}

void DarkQCD::InitDarkQCDVertices() {
  Settings& s = Settings::GetMainSettings();
  if (!Flavour(kf_dark_g).IsOn()) return;
  m_dec_dark_g4 = false; //s["DECOMPOSE_4G_VERTEX"].Get<int>();
  Kabbala g3("g_3_dark",sqrt(4.*M_PI*ScalarConstant("alpha_Dark")));
  Kabbala cpl0=g3*Kabbala("i",Complex(0.,1.));
  Flavour flav;
  for (long int i=4900101;i<4900103;++i) {
    Flavour flav((kf_code)i);
    if (flav.IsOn()) {
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle(flav.Bar());
      m_v.back().AddParticle(flav);
      m_v.back().AddParticle(Flavour(kf_dark_g));
      m_v.back().Color.push_back(Color_Function(cf::T,3,2,1));
      m_v.back().Lorentz.push_back("FFV");
      m_v.back().cpl.push_back(cpl0);
      /* m_v.back().order[0]=1; */
      m_v.back().order.push_back(1);
    }
  }
  // flav = Flavour(kf_Dv);
  // if (flav.IsOn()) {
  //   m_v.push_back(Single_Vertex());
  //   m_v.back().AddParticle(flav.Bar());
  //   m_v.back().AddParticle(flav);
  //   m_v.back().AddParticle(Flavour(kf_gv));
  //   m_v.back().Color.push_back(Color_Function(cf::T,3,2,1));
  //   m_v.back().Lorentz.push_back("FFV");
  //   m_v.back().cpl.push_back(cpl0);
  //   /* m_v.back().order[0]=1; */
  //   m_v.back().order.push_back(1);
  // }
  Kabbala cpl1=-g3;
  m_v.push_back(Single_Vertex());
  for (size_t i(0);i<3;++i) m_v.back().AddParticle(Flavour(kf_dark_g));
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
  for (size_t i(0);i<4;++i) m_v.back().AddParticle(Flavour(kf_dark_g));
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

void DarkQCD::AddHadronDecays(PHASIC::Decay_Map * dmap) {
  msg_Out()<<METHOD<<" for "<<m_hadrons.size()<<" hadron types.\n";
  p_ms  = dmap->GetMassSelector();
  p_mix = (dynamic_cast<HADRONS::Hadron_Decay_Map *>(dmap))->MixingHandler();
  for (set<Flavour,SortByMass>::iterator hit=m_hadrons.begin();hit!=m_hadrons.end();hit++) {
    if ((*hit).Stable()==1) {
      msg_Out()<<"  * "<<(*hit)<<" is stable.  Continue.\n"; continue;
    }
    else {
      msg_Out()<<"  * "<<(*hit)<<" decays.  Construct decay table.\n";
      if ((*hit)==Flavour(kf_dark_pi11) || (*hit)==Flavour(kf_dark_pi21) ||
	  (*hit)==Flavour(kf_dark_eta22)) {
	HADRONS::Hadron_Decay_Table * hdt = DecaysForPis((*hit));
	if (hdt) (*dmap)[(*hit)] = hdt;
      }
      else if ((*hit)==Flavour(kf_dark_rho11) || (*hit)==Flavour(kf_dark_rho21) ||
	       (*hit)==Flavour(kf_dark_omega22)) {
	HADRONS::Hadron_Decay_Table * hdt = DecaysForRhos((*hit));
	if (hdt) (*dmap)[(*hit)] = hdt;
      }
    }
  }
}

HADRONS::Hadron_Decay_Table * DarkQCD::DecaysForRhos(const Flavour & flav) {
  HADRONS::Hadron_Decay_Table * hdt = nullptr;
  double Gamma = 0.;
  double M     = flav.Mass(), M2 = sqr(M);
  std::vector<Flavour> pions;
  if (flav==Flavour(kf_dark_rho11) || flav==Flavour(kf_dark_omega22)) {
    double m   = Flavour(kf_dark_pi21).Mass();
    if (M>2.*m) {
      double m2    = sqr(m);
      double Ampl2 = sqr(flav==Flavour(kf_dark_rho11) ? m_gRhoPiPi : m_gOmegaPiPi)/3. * (M2-4.*m2);
      Gamma = GammaTwoBody(Ampl2,M2,m2,m2);
      pions.push_back(Flavour(kf_dark_pi21));
      pions.push_back(Flavour(kf_dark_pi21).Bar());
    }
    else {
      msg_Error()<<"Error in "<<METHOD<<":\n"
		 <<"   Mass of dark rho below 2-pion threshold.\n"
		 <<"   Todo: Will assume decay into fermions as 'mixing' with the Z'.\n"
		 <<"   Also: will assume the mixing parameter is given by g^2.\n";
      exit(1);
    }
  }
  if (flav==Flavour(kf_dark_rho21)) {
    double m2 = Flavour(kf_dark_pi21).Mass();
    double m1 = Flavour(kf_dark_pi11).Mass();
    if (M>m1+m2) {
      double m12   = sqr(m1);
      double m22   = sqr(m2);
      double Ampl2 = sqr(m_gRhoPiPi)/3. * (sqr(M2-m12-m22)-4.*m12*m22)/M2;
      Gamma = GammaTwoBody(Ampl2,M2,m12,m22);
      pions.push_back(Flavour(kf_dark_pi21));
      pions.push_back(Flavour(kf_dark_pi11));
    }
    else {
      msg_Error()<<"Error in "<<METHOD<<":\n"
		 <<"   Mass of dark rho below 2-pion threshold.\n"
		 <<"   Todo: Will assume decay into fermions as 'mixing' with the Z'.\n"
		 <<"   Also: will assume the mixing parameter is given by g^2.\n";
      exit(1);
    }
  }
  if (Gamma>0.) {
    const Scoped_Settings &s = Settings::GetMainSettings()["Dark_QCD"];
    HADRONS::GeneralModel md;
    flav.SetWidth(Gamma);
    hdt = new HADRONS::Hadron_Decay_Table(flav,p_ms,p_mix);
    HADRONS::Hadron_Decay_Channel * hdc = new HADRONS::Hadron_Decay_Channel(flav,p_ms);
    while (!pions.empty()) {
      hdc->AddDecayProduct(pions.back());
      pions.pop_back();
    }
    hdc->SetWidth(Gamma);
    hdc->SetDeltaWidth(0.);
    hdc->Initialise(s,md);
    hdt->AddDecayChannel(hdc);
  }
  msg_Out()<<(*hdt)<<"\n";
  return hdt;
}

HADRONS::Hadron_Decay_Table * DarkQCD::DecaysForPis(const Flavour & pion) {
  HADRONS::Hadron_Decay_Table * hdt = nullptr;
  std::map<ATOOLS::Flavour,double> gammas;
  double Gamma = 0.;
  double mpi2  = sqr(pion.Mass());
  if (pion==Flavour(kf_dark_pi11) || pion==Flavour(kf_dark_eta22)) {
    Flavour Zp = Flavour(kf_Z0_2);
    if (Zp.IsOn()) {
      double MZ2   = sqr(Zp.Mass());
      double pref  = pion==Flavour(kf_dark_pi11) ? sqr(m_fPi/MZ2) : sqr(m_fEta/MZ2);
      msg_Out()<<METHOD<<" for "<<pion<<", mpi2 = "<<mpi2<<", mZ2 = "<<MZ2<<", fpi = "<<m_fPi<<".\n";
      // assume <0|J^mu|pi> = f_pi * P^mu
      // assume Z' being off-shell and (we are only interested in
      // approximate gamma's to resolve the correct branching ratios.
      // ampl = (f_pi * P^mu) * (-g_{mu nu}+P_mu P_nu/MZ^2)/(P^2-MZ^2) *
      //        bar{u}(p_1) gamma^nu p_2
      //      = - [f_pi (m1+m2) (1-mpi^2/MZ^2)] / (mpi^2-MZ^2)
      for (short int i=1;i<=16;++i) {
	if (i==7) i = 11;
	Flavour flav((kf_code)i);
	if (flav.IsOn()) {
	  double m1 = flav.HadMass(), m2 = m1;
	  double g  = 0., ampl = 0.;
	  if (sqr(m1+m2)>mpi2) continue;
	  for (size_t j=0;j<m_v.size();j++) {
	    Single_Vertex v = m_v[j];
	    if (v.in[0]==flav.Bar() && v.in[1]==flav && v.in[2]==Zp) {
	      g = abs(v.cpl[0].Value());
	      if (flav.IsQuark()) g *= 3.;
	      double Ampl2 = pref*sqr(g)*sqr(m1+m2);
	      Gamma += gammas[flav] = GammaTwoBody(Ampl2,mpi2,sqr(m1),sqr(m2));
	      msg_Out()<<"Gotcha. "<<flav<<"("<<flav.IsQuark()<<"): m1 = "<<m1<<", g = "<<g<<" (ampl = "<<Ampl2<<") "
		       <<"--> gamma = "<<gammas[flav]<<" --> Gamma = "<<Gamma<<".\n";
	    }
	  }
	}
      }
    }
  }
  if (Gamma>0.) {
    const Scoped_Settings &s = Settings::GetMainSettings()["Dark_QCD"];
    HADRONS::GeneralModel md;
    pion.SetWidth(Gamma);
    hdt = new HADRONS::Hadron_Decay_Table(pion,p_ms,p_mix);
    for (std::map<ATOOLS::Flavour,double>::iterator cit=gammas.begin();cit!=gammas.end();cit++) {
      HADRONS::Hadron_Decay_Channel * hdc = new HADRONS::Hadron_Decay_Channel(pion,p_ms);
      hdc->AddDecayProduct(cit->first);
      hdc->AddDecayProduct(cit->first.Bar());
      hdc->SetWidth(cit->second);
      hdc->SetDeltaWidth(0.);
      hdc->Initialise(s,md);
      hdt->AddDecayChannel(hdc);
    }
  }
  msg_Out()<<(*hdt)<<"\n";
  return hdt;
}

double DarkQCD::GammaTwoBody(const double & Ampl2,const double & M2,
			     const double & m12,const double & m22) {
  return Ampl2/(16.*M_PI*sqrt(M2))*sqrt(sqr(M2-(m12+m22))-4.*m12*m22)/M2;
}
