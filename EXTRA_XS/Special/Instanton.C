#include "PHASIC++/Process/External_ME_Args.H"
#include "PHASIC++/Channels/Rambo.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include "MODEL/Main/Model_Base.H"
#include "MODEL/Main/Running_AlphaS.H"

#include "EXTRA_XS/Main/ME2_Base.H"

using namespace EXTRAXS;
using namespace MODEL;
using namespace ATOOLS;
using namespace PHASIC;
using namespace std;


/*
   In all the differential cross sections the factor 1/16 Pi is cancelled
   by the factor 4 Pi for each alpha. Hence one Pi remains in the game.
*/

namespace EXTRAXS {
  struct xsec_data {
    double m_Ehat, m_rho, m_Ngluons, m_sigmahat;
    xsec_data(const double & Ehat, const double & rho, const double & Ngluons,
	      const double & sigmahat) :
      m_Ehat(Ehat), m_rho(rho), m_Ngluons(Ngluons), m_sigmahat(sigmahat) {}
  };

  class Data_Table {
  private:
    std::map<double, xsec_data *> m_data;
    double      m_Ehatmin, m_Ehatmax;
    double      m_rho, m_Ngluons, m_sigmahat;
    xsec_data * m_lower, * m_upper;
    bool ReadTable();
    void ConstructDefaultTable();
  public:
    Data_Table();
    ~Data_Table();
    void Output();
    bool Interpolate(const double & E);
    void Test(const double & E);
    const double & Ehatmin()  const { return m_Ehatmin; }
    const double & Ehatmax()  const { return m_Ehatmax; }
    const double & Rho()      const { return m_rho; }
    const double & Ngluons()  const { return m_Ngluons; }
    const double & Sigmahat() const { return m_sigmahat; }
    const double & E(int & i) {
      std::map<double, xsec_data*>::iterator dit = m_data.begin();
      while ((i--)>0) dit++;
      return dit->first;
    }
  };

  Data_Table::Data_Table() {
    if (!ReadTable()) ConstructDefaultTable();
    m_Ehatmin = m_data.begin()->first;
    m_Ehatmax = m_data.rbegin()->first;
    for (std::map<double, xsec_data*>::iterator dit=m_data.begin();
	 dit!=m_data.end();dit++)
      dit->second->m_sigmahat = dit->second->m_sigmahat/rpa->Picobarn();
  }

  Data_Table::~Data_Table() {
    while (!m_data.empty()) {
      delete m_data.begin()->second;
      m_data.erase(m_data.begin());
    }
    m_data.clear();
  }

  bool Data_Table::ReadTable() {
    Settings& s = Settings::GetMainSettings();
    s.DeclareMatrixSettingsWithEmptyDefault({ "INSTANTON_XSECS" });
    const auto allxs = s["INSTANTON_XSECS"].GetMatrix<double>();
    if (allxs.empty()) {
      msg_Info()<<"Warning in "<<METHOD<<":\n"
		<<"   did not find matrix of cross sections under tag INSTANTON_XSECS\n"
		<<"   will read defaults from 1911.09726v1.\n";
      return false;
    }
    for (const auto xs : allxs) {
      /*
	read data from file in ASCII format: in each line we have
	\sqrt{s} (GeV) , 1/\rho (GeV) , \alpha_s(1/rho) ,  number of gluons ,
	Cross-section (pb)
      */
      if (xs.size()!=5) {
	msg_Error()<<"Error in "<<METHOD<<":\n"
		   <<"   Badly formatted xs entry in tag INSTANTON_XSECS, \n"
		   <<"   needs 5 doubles (E, 1/rho, alphaS, Ngluons, xsec).\n"
		   <<"   will read defaults from 1911.09726v1.\n";
	return false;
      }
      msg_Tracking()<<"   xs("<<xs[0]<<") = "<<xs[4]<<" pb.\n";
      m_data[xs[0]] = new xsec_data{ xs[0], xs[1], xs[3], xs[4] };
    }
    return true;
  }

  void Data_Table::ConstructDefaultTable() {
  /*
    {\sqrt{s} (GeV) , 1 / \rho (GeV) , \alpha_s (at 1/rho) ,  number of gluons ,
    Cross-section (pb)} where the cross-section allows for N_f=4 plus N_f=5 if
    it’s above the threshold (so always 8 to 10 fermions in addition to n_g gluons
    in the final state).

    You can see that the largest cross-section is 14.5 milli barn at 10 GeV
    (with alpha_s = 0.415 still kind of reasonable) and dropping to 4 pico barns
    at 40 GeV

    This is the original data from the first version of 1911.09726

    {{10.6853, 0.989378, 0.415544, 4.58901, 1.45813*10^10},
    {11.3923, 1.03566, 0.40533, 4.67934, 1.05266*10^10},
    {13.3679, 1.16243,  0.38164, 4.90476, 4.51405*10^9},
    {15.6816, 1.3068, 0.360291, 5.12963, 1.85274*10^9},
    {22.9076, 1.76212, 0.315257, 5.43739, 1.76977*10^8},
    {29.6526, 2.11804, 0.292739, 6.02483, 3.55261*10^7},
    {40.769, 2.71793, 0.266897, 6.47057, 3.99487*10^6}}
    {56.0679, 3.50425, 0.244871, 6.91549, 397757.},
    {61.8378, 3.63752, 0.223469, 7.28038, 113207.},
    {89.628, 4.97933, 0.205845, 7.67006, 3876.77},
    {118.028, 6.212, 0.195011, 8.24762, 333.886},
    {174.399, 8.71997, 0.180449, 8.60356, 8.68739},
    {246.887, 11.7565, 0.169311, 9.04486, 0.348676},
    {349.853, 15.9024, 0.159366, 9.48648, 0.0140647},
    {496.283, 21.5775, 0.150438, 9.9284, 0.000571738}}
    {704.764, 29.3652, 0.142384, 10.3706, 0.0000232145},
    {1001.82, 40.0727, 0.135088, 10.8133, 9.29353*10^-7},
    {1425.6, 54.8307, 0.128449, 11.2564, 3.61946*10^-8},
    {2030.63, 75.2085, 0.122387, 11.7001, 1.36042*10^-9},
    {2895.48, 103.41, 0.116831, 12.1441, 4.83086*10^-11},
    {4132.81, 142.511, 0.111723, 12.5887, 1.68134*10^-12}}
  */
    m_data[10.6853] = new xsec_data(10.6853, 0.989378, 4.58901, 1.45813e10);
    m_data[11.3932] = new xsec_data(11.3923, 1.035660, 4.67934, 1.05266e10);
    m_data[13.3679] = new xsec_data(13.3679, 1.162430, 4.90476, 4.51405e9);
    m_data[15.6816] = new xsec_data(15.6816, 1.306800, 5.12963, 1.85274e9);
    m_data[22.9076] = new xsec_data(22.9076, 1.762120, 5.43739, 1.76977e8);
    m_data[29.6526] = new xsec_data(29.6526, 2.118040, 6.02483, 3.55261e7);
    m_data[40.7690] = new xsec_data(40.7690, 2.717930, 6.47057, 3.99487e6);
    m_data[56.0679] = new xsec_data(56.0679, 3.504250, 6.91549, 397757.);
    m_data[61.8378] = new xsec_data(61.8378, 3.637520, 7.28038, 113207.);
    m_data[89.6280] = new xsec_data(89.6280, 4.979330, 7.67006, 3876.77);
    m_data[118.028] = new xsec_data(118.028, 6.212000, 8.24762, 333.886);
    m_data[174.399] = new xsec_data(174.399, 8.719970, 8.60356, 8.68739);
    m_data[246.887] = new xsec_data(246.887, 11.75650, 9.04486, 0.348676);
    m_data[349.853] = new xsec_data(349.853, 15.90240, 9.48648, 0.0140647);
    m_data[496.283] = new xsec_data(496.283, 21.57750, 9.92840, 0.000571738);
    m_data[704.764] = new xsec_data(704.764, 29.36520, 10.3706, 0.0000232145);
    m_data[1001.82] = new xsec_data(1001.82, 40.07270, 10.8133, 9.29353e-7);
    m_data[1425.60] = new xsec_data(1425.60, 54.83070, 11.2564, 3.61946e-8);
    m_data[2030.63] = new xsec_data(2030.63, 75.20850, 11.7001, 1.36042e-9);
    m_data[2895.48] = new xsec_data(2895.48, 103.4100, 12.1441, 4.83086e-11);
    m_data[4132.81] = new xsec_data(4132.81, 142.5110, 12.5887, 1.68134e-12);
  }

  bool Data_Table::Interpolate(const double & E) {
    m_rho = m_sigmahat = m_Ngluons = 0.;
    if (E>m_Ehatmax || E<m_Ehatmin) {
      msg_Debugging()<<METHOD<<" yields false for E = "<<E<<".\n";
      return false;
    }
    std::map<double, xsec_data*>::iterator dit;
    for (dit=m_data.begin();dit!=m_data.end();dit++) if (dit->first>E) break;
    m_upper = dit->second;
    dit--;
    m_lower = dit->second;
    double lx = (m_upper->m_Ehat-E)/(m_upper->m_Ehat-m_lower->m_Ehat);
    double ux = (m_lower->m_Ehat-E)/(m_lower->m_Ehat-m_upper->m_Ehat);
    m_rho           = (m_upper->m_rho*ux      + m_lower->m_rho*lx);
    m_Ngluons       = (m_upper->m_Ngluons*ux  + m_lower->m_Ngluons*lx);
    m_sigmahat      = (m_upper->m_sigmahat*ux + m_lower->m_sigmahat*lx);
    return true;
  }

  void Data_Table::Output() {
    msg_Out()<<"Instanton partonic cross sections:\n"
	     <<"   with sqrt(s') in ["<<m_Ehatmin<<", "<<m_Ehatmax<<"]\n"
	     <<"--------------------------------------------------\n"
	     <<"E'[GeV]:   1/rho[GeV]    <Ngluons>   sigma\n";
    for (std::map<double, xsec_data*>::iterator dit=m_data.begin();
	 dit!=m_data.end();dit++)
      msg_Out()<<dit->first<<": "<<dit->second->m_rho<<" "
	       <<dit->second->m_Ngluons<<" "
	       <<(dit->second->m_sigmahat*rpa->Picobarn())<<"\n";
    msg_Out()<<"--------------------------------------------------\n";
  }

  void Data_Table::Test(const double & E) {
    Output();
    Interpolate(E);
    msg_Out()<<"For E = "<<E<<" GeV: sigma' = "<<Sigmahat()<<" 1/GeV^2 = "
	     <<Sigmahat()*rpa->Picobarn()<<" pb, "
	     <<"1/rho = "<<Rho()<<" GeV, "
	   <<"<Ngluons> = "<<Ngluons()<<".\n";
  }

  struct instantonScale {
    enum code {
      rho           = 1,
      Ehat          = 2,
      Ehat_by_sqrtN = 3
    };
  };

  class XS_instanton : public ME2_Base {
  private:
    Data_Table  m_data;
    double      m_Ehatmin, m_Ehatmax, m_norm, m_S, m_Ehat;
    double      m_tthreshold, m_bthreshold, m_cthreshold;
    double      m_Ngluons_factor, m_sigmahat_factor, m_alphaS_factor, m_scale_factor;
    instantonScale::code m_scalechoice;
    double      m_Ecms, m_threshold, m_mean_Ngluons;
    size_t      m_nquarks, m_ngluons, m_includeQ;
    MODEL::Running_AlphaS * p_alphaS;
    std::vector<double>     m_masses;
    void   Initialise();
    bool   DefineFlavours();
    double FixScale() const;
    double AlphaSModification();
    size_t NumberOfGluons();
    bool   DistributeMomenta();
    bool   MakeColours();
    void   Test();
  public:
    XS_instanton(const External_ME_Args& args);
    ~XS_instanton() {} // if (p_data) { delete p_data; p_data = NULL; }}
    double operator()(const ATOOLS::Vec4D_Vector& mom) override;
    double RelativeRenormalisationScaleWeight(double scalefactor) const override;
    bool   SetColours(const Vec4D_Vector& mom) override;
    bool   FillFinalState(const std::vector<ATOOLS::Vec4D> & mom) override;

    // Report that this class has a non-standard AlphaS dependency, but offers
    // RelativeRenormalisationScaleWeight to calculate it.
    int OrderQCD(const int&) const override { return 98; };
  };
}

XS_instanton::XS_instanton(const External_ME_Args& args)
  : ME2_Base(args),
    m_S(sqr(rpa->gen.Ecms())),
    m_norm(1./36.),
    m_bthreshold(100.), m_cthreshold(20.),
    m_Ngluons_factor(1.), m_sigmahat_factor(1.),
    m_scalechoice(instantonScale::rho)
{
  DEBUG_INFO("now entered EXTRAXS::XS_instanton ...");
  Settings& s       = Settings::GetMainSettings();
  double Ehatmin    = Max(1.,m_data.Ehatmin());
  double Ehatmax    = Min(rpa->gen.Ecms(),m_data.Ehatmax());
  m_Ehatmin         = s["INSTANTON_MIN_MASS"].SetDefault(20.).Get<double>();
  m_Ehatmax         = s["INSTANTON_MAX_MASS"].SetDefault(rpa->gen.Ecms()).Get<double>();
  m_sprimemin       = sqr(m_Ehatmin);
  m_sprimemax       = sqr(m_Ehatmax);
  if (m_Ehatmin<Ehatmin || m_Ehatmax>Ehatmax) {
    msg_Error()<<"WARNING in "<<METHOD<<":\n"
	       <<"   mass range of simulation not fully captured by data:\n";
    if (m_Ehatmin<Ehatmin) {
      msg_Error()<<"   demand minimal instanton mass below smallest energy in data:\n"
		 <<"   "<<m_Ehatmin<<" < "<<Ehatmin
		 <<" -- this could be a problem due to critical extrapolation.\n";
    }
    if (m_Ehatmax>Ehatmax) {
      msg_Error()<<"   demand maximal instanton mass above largest energy in data:\n"
		 <<"   "<<m_Ehatmax<<" > "<<Ehatmax
		 <<" -- this should not become a problem.\n";
    }
  }
  m_Ngluons_factor  = s["INSTANTON_NGLUONS_MODIFIER"].SetDefault(1.0).Get<double>();
  m_tthreshold      = (s["INSTANTON_T_PRODUCTION_THRESHOLD"].SetDefault(1000.0).
			Get<double>());
  m_bthreshold      = (s["INSTANTON_B_PRODUCTION_THRESHOLD"].SetDefault(100.0).
			Get<double>());
  m_cthreshold      = (s["INSTANTON_C_PRODUCTION_THRESHOLD"].SetDefault(20.0).
			Get<double>());
  m_includeQ        = s["INSTANTON_INCLUDE_QUARKS"].SetDefault(5).Get<int>();
  m_sigmahat_factor = s["INSTANTON_SIGMAHAT_MODIFIER"].SetDefault(1.0).Get<double>();
  m_alphaS_factor   = s["INSTANTON_ALPHAS_FACTOR"].SetDefault(1.0).Get<double>();
  string choice     = s["INSTANTON_SCALE_CHOICE"].SetDefault("1/rho").Get<string>();
  if      (choice=="shat")   m_scalechoice = instantonScale::Ehat;
  else if (choice=="shat/N") m_scalechoice = instantonScale::Ehat_by_sqrtN;
  m_scale_factor    = s["INSTANTON_SCALE_FACTOR"].SetDefault(1.).Get<double>();


  p_alphaS = dynamic_cast<Running_AlphaS *>(s_model->GetScalarFunction("alpha_S"));
  m_hasinternalscale = true;
  msg_Tracking()<<METHOD<<" for instanton production in the energy range "
		<<"["<<m_Ehatmin<<", "<<m_Ehatmax<<"]\n"
		<<"   Ngluons factor = "<<m_Ngluons_factor<<", "
		<<"sigmahat factor = "<<m_sigmahat_factor<<".\n";
}

double XS_instanton::operator()(const Vec4D_Vector& momenta) {
  double shat = momenta[2].Abs2();
  m_Ehat = sqrt(shat);
  if (m_Ehat<m_Ehatmin || m_Ehat>m_Ehatmax ||
      !m_data.Interpolate(m_Ehat)) return 0.;
  m_internalscale = Max(FixScale(), 2.);
  // have to multiply with the norm and the inverse external flux
  double xsec = m_sigmahat_factor * m_data.Sigmahat() * (2.*shat) * m_norm;
  return AlphaSModification() * xsec;
}

double XS_instanton::FixScale() const {
  return m_scale_factor *
    (m_scalechoice==instantonScale::Ehat?m_Ehat:
     (m_scalechoice==instantonScale::Ehat_by_sqrtN)?m_Ehat/sqrt(m_data.Ngluons()):
     m_data.Rho());
}

double XS_instanton::AlphaSModification() {
  if (dabs(m_alphaS_factor-1.)<1.e-3) return 1.;
  return pow(m_alphaS_factor,2.*p_alphaS->Beta0(sqr(m_data.Rho())));
}

double XS_instanton::RelativeRenormalisationScaleWeight(double scalefactor) const {
  // This is based on the dominant scale dependence of the instanton cross section,
  // cf. Eq. (2.30) of 1911.09726
  if (dabs(scalefactor-1.)<1.e-3) return 1.;
  double scale2 = sqr(FixScale());
  return ( exp(-4.*M_PI * (1./(*p_alphaS)(scalefactor*scale2)-1./(*p_alphaS)(scale2))) *
	   pow(sqrt(scalefactor)/m_alphaS_factor,2.*4.*p_alphaS->Beta0(scalefactor*scale2)) );
}


bool XS_instanton::FillFinalState(const std::vector<Vec4D> & mom) {
  m_Ehat = sqrt(mom[2].Abs2());
  if (m_Ehat<m_Ehatmin || m_Ehat>m_Ehatmax ||
      !m_data.Interpolate(m_Ehat)) return false;
  Poincare boost(mom[2]);
  m_internalscale = Max(FixScale(), 2.);
  //msg_Out()<<METHOD<<": scale = "<<m_internalscale<<" -> "<<sqr(m_internalscale)<<"\n";
  size_t trials=1000;
  while ((trials--)>0) {
    m_mean_Ngluons = m_data.Ngluons();
    Vec4D sum      = -mom[2];
    if (DefineFlavours() && DistributeMomenta() && MakeColours()) {
      for (size_t i=0;i<m_flavours.size();i++) {
	boost.BoostBack(m_momenta[i]);
	sum += m_momenta[i];
      }
      if (dabs((mom[2]-sum).Abs2())<1.e-6 && dabs(mom[2][0]-sum[0])<1.e-6) {
	//msg_Out()<<METHOD<<" for "<<m_flavours.size()<<"\n";
	return true;
      }
    }
  }
  return false;
}

bool XS_instanton::DefineFlavours() {
  m_flavours.clear();
  double totmass = 0.;
  for (size_t i=0;i<2;i++) {
    totmass += m_flavs[i].Mass(true);
    m_flavours.push_back(m_flavs[i]);
  }
  m_masses.clear();
  m_nquarks = 0;
  for (size_t i=1;i<6;i++) {
    Flavour flav = Flavour(i);
    if (i>m_includeQ)                               continue;
    if (flav==Flavour(kf_b) && m_Ehat<m_bthreshold) continue;
    if (flav==Flavour(kf_c) && m_Ehat<m_cthreshold) continue;
    if (totmass+2.*flav.Mass(true)>m_Ehat) break;
    if (flav.Bar()!=m_flavs[0] && flav.Bar()!=m_flavs[1]) {
      m_flavours.push_back(flav);
      m_nquarks++;
      totmass += flav.Mass(true);
    }
    flav = flav.Bar();
    if (flav.Bar()!=m_flavs[0] && flav.Bar()!=m_flavs[1]) {
      m_flavours.push_back(flav);
      m_nquarks++;
      totmass += flav.Mass(true);
    }
  }
  do { m_ngluons = NumberOfGluons(); } while (m_ngluons>=30-m_nquarks);
  Flavour flav   = Flavour(kf_gluon);
  for (size_t i=0;i<m_ngluons;i++)  m_flavours.push_back(flav);
  return true;
}

size_t XS_instanton::NumberOfGluons() {
  return ran->Poissonian(m_Ngluons_factor * m_mean_Ngluons);
}

bool XS_instanton::DistributeMomenta() {
  m_momenta.clear();
  double totmass = 0., mass;
  for (size_t i=0;i<m_flavours.size();i++) {
    totmass += mass = m_flavours[i].Mass(true);
    m_masses.push_back(mass);
  }
  if (totmass>m_Ehat) {
    msg_Error()<<"Error in "<<METHOD<<" not enough energy (Ecms = "<<m_Ehat<<") "
	       <<" in instanton to produce "<<(m_flavours.size()-2)<<" "
	       <<"partons in its decay,\n"
	       <<"   total partonic FS mass is "<<mass<<"\n";
    return false;
  }
  Rambo rambo(2,m_masses);
  m_momenta = rambo.GeneratePoint(m_Ehat);
  return true;
}

bool XS_instanton::MakeColours() {
  for (size_t i=0;i<m_colours.size();i++) m_colours[i].clear(); m_colours.clear();
  m_colours.resize(m_flavours.size());
  vector<size_t> cols[2];
  for (size_t i=0;i<m_flavours.size();i++) {
    Flavour flav = m_flavours[i];
    if ((flav.IsQuark() && !flav.IsAnti()) || flav.IsGluon())
      cols[i<2?1:0].push_back(i);
    if ((flav.IsQuark() &&  flav.IsAnti()) || flav.IsGluon())
      cols[i<2?0:1].push_back(i);
    m_colours[i].resize(2);
  }
  size_t pos[2], parts[2], colindex = 500;
  do {
    do {
      for (size_t i=0;i<2;i++) {
	pos[i]   = size_t(0.999999999*cols[i].size()*ran->Get());
	parts[i] = cols[i][pos[i]];
      }
    } while (parts[0]==parts[1]);
    m_colours[parts[0]][0] = m_colours[parts[1]][1] = colindex++;
    for (size_t i=0;i<2;i++) { cols[i].erase(find(cols[i].begin(),cols[i].end(),parts[i])); }
  } while (cols[0].size()>1);
  if (cols[0][0]!=cols[1][0]) {
    m_colours[cols[0][0]][0] = m_colours[cols[1][0]][1] = colindex++;
  }
  else {
    if (m_colours[2][0]!=0) {
      m_colours[cols[0][0]][0] = m_colours[2][0];
      m_colours[2][0] = m_colours[cols[1][0]][1] = colindex++;
    }
    else if (m_colours[2][1]!=0) {
      m_colours[cols[1][0]][1] = m_colours[2][1];
      m_colours[2][1] = m_colours[cols[0][0]][0] = colindex++;
    }
  }
  for (size_t i=0;i<2;i++) {
    size_t help     = m_colours[i][0];
    m_colours[i][0] = m_colours[i][1];
    m_colours[i][1] = help;
  }
  return true;
}

void XS_instanton::Test() {
  for (size_t i=0;i<5;i++) {
    int step = 2*i;
    double E = m_data.E(step);
    m_data.Interpolate(E);
    m_mean_Ngluons = m_data.Ngluons();
    long int maxtrials = 1000000, ngluons = 0, totn = 0;
    for (long int trials=0;trials<maxtrials;trials++) {
      totn += ngluons = NumberOfGluons();
    }
    msg_Out()<<"Run for E = "<<E<<": <ngluons> = "<<m_mean_Ngluons
	     <<" --> "<<double(totn)/double(maxtrials)<<" generated.\n";
  }
}

bool XS_instanton::SetColours(const Vec4D_Vector& mom) {
  THROW(fatal_error,"XS_instanton::SetColours should never be called.");
}

DECLARE_TREEME2_GETTER(EXTRAXS::XS_instanton,"XS_instanton")
Tree_ME2_Base *ATOOLS::Getter<PHASIC::Tree_ME2_Base,PHASIC::External_ME_Args,EXTRAXS::XS_instanton>::
operator()(const External_ME_Args& args) const
{
  if (MODEL::s_model->Name()!="SM") return NULL;
  const Flavour_Vector fl=args.Flavours();
  if (fl.size()!=3) return NULL;
  if (!(fl[0].Strong() && fl[1].Strong() &&
	fl[2]==Flavour(kf_instanton))) return NULL;
  return new XS_instanton(args);
}
