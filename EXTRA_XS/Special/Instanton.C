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
  public:
    Data_Table();
    ~Data_Table();
    void Output();
    bool Interpolate(const double & E);
    void Test();
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

  /*
    {\sqrt{s} (GeV) , 1 / \rho (GeV) , \alpha_s (at 1/rho) ,  number of gluons , 
    Cross-section (pb)} where the cross-section allows for N_f=4 plus N_f=5 if it’s above 
    the threshold (so always 8 to 10 fermions in addition to n_g gluons in the final state)
    
    You can see that the largest cross-section is 14.5 milli barn at 10 GeV 
    (with alpha_s = 0.415 still kind of reasonable) and dropping to 4 pico barns at 40 GeV
  */
  Data_Table::Data_Table() {
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

  void Data_Table::Test() {
    Output();
    Interpolate(100.);
    msg_Out()<<"For E = 100 GeV: sigma' = "<<Sigmahat()<<" 1/GeV^2 = "
	     <<Sigmahat()*rpa->Picobarn()<<" pb, "
	     <<"1/rho = "<<Rho()<<" GeV, "
	   <<"<Ngluons> = "<<Ngluons()<<".\n";
    exit(1);
  }
  
  class XS_instanton : public ME2_Base {  // == XS_ffbar_ee but not XS_ffbar_f'fbar' !
  private:
    Data_Table              m_data;
    double                  m_Ehatmin, m_Ehatmax, m_norm, m_S, m_Ehat;
    double                  m_Ngluons_factor, m_sigmahat_factor;
    std::string             m_scalechoice, m_Qthreshold;
    double                  m_Ecms, m_threshold, m_mean_Ngluons;
    size_t                  m_nquarks, m_ngluons;
    MODEL::Running_AlphaS * p_alphaS;
    std::vector<double>     m_masses;
    bool   DefineFlavours();
    size_t NumberOfGluons();
    bool   DistributeMomenta();
    bool   MakeColours();
    void   Test();
  public:
    XS_instanton(const External_ME_Args& args);
    ~XS_instanton() {}
    double operator()(const ATOOLS::Vec4D_Vector& mom);
    bool   SetColours(const Vec4D_Vector& mom);
    bool   FillFinalState(const std::vector<ATOOLS::Vec4D> & mom);
  };
}

XS_instanton::XS_instanton(const External_ME_Args& args)
  : ME2_Base(args),
    m_Ehatmin(Max(1.,m_data.Ehatmin())),
    m_Ehatmax(Min(rpa->gen.Ecms(),m_data.Ehatmax())),
    m_S(sqr(rpa->gen.Ecms())),
    m_norm(1./36.),
    m_Ngluons_factor(1.), m_sigmahat_factor(1.)
{
  p_alphaS    = dynamic_cast<Running_AlphaS *>(s_model->GetScalarFunction("alpha_S"));
  Settings& s = Settings::GetMainSettings();
  DEBUG_INFO("now entered EXTRAXS::XS_instanton ...");
  m_Ngluons_factor   = s["INSTANTON_NGLUONS_MODIFIER"].SetDefault(1.0).Get<double>();
  m_sigmahat_factor  = s["INSTANTON_SIGMAHAT_MODIFIER"].SetDefault(1.0).Get<double>();
  m_Ehatmin          = Max(s["INSTANTON_MIN_MASS"].SetDefault(m_Ehatmin).Get<double>(),
			   m_Ehatmin);
  m_Ehatmax          = Min(s["INSTANTON_MAX_MASS"].SetDefault(m_Ehatmax).Get<double>(),
			   m_Ehatmax);
  m_sprimemin        = sqr(m_Ehatmin);
  m_sprimemax        = sqr(m_Ehatmax);
  m_scalechoice      = s["INSTANTON_SCALE_CHOICE"].SetDefault("1/rho").Get<string>();
  m_Qthreshold       = s["INSTANTON_QUARK_THRESHOLD"].SetDefault("shat").Get<string>();
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
  m_internalscale = (m_scalechoice==std::string("shat"))?m_Ehat:m_data.Rho();
  //if (m_internalscale<2.) m_internalscale = 2.;
  // have to multiply with the norm and the inverse external flux
  double xsec = m_sigmahat_factor * m_data.Sigmahat() * (2.*shat) * m_norm;
  return xsec;
}

bool XS_instanton::FillFinalState(const std::vector<Vec4D> & mom) {
  m_Ehat = sqrt(mom[2].Abs2());
  if (m_Ehat<m_Ehatmin || m_Ehat>m_Ehatmax ||
      !m_data.Interpolate(m_Ehat)) return false;
  m_internalscale = (m_scalechoice==std::string("shat"))?m_Ehat:m_data.Rho();
  m_threshold     = (m_Qthreshold==std::string("shat"))?m_Ehat:m_data.Rho();
  m_mean_Ngluons  = m_data.Ngluons();
  Poincare boost(mom[2]);
  
  if (DefineFlavours() && DistributeMomenta() && MakeColours()) {
    for (size_t i=0;i<m_flavours.size();i++) boost.BoostBack(m_momenta[i]);
    //Vec4D check = mom[2];
    for (size_t i=0;i<m_flavours.size();i++) {
      boost.BoostBack(m_momenta[i]);
      //if (i>=2) check -= m_momenta[i];
      //msg_Out()<<"   * "<<i<<": "<<m_flavours[i]<<" "
      //	       <<"["<<p_colours[i][0]<<" "<<p_colours[i][1]<<"] "
      //       <<m_momenta[i]<<"\n";
    }
    //msg_Out()<<"   *** 4 mom: "<<mom[2]<<" --> "<<check<<".\n";
  }
  return true;
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
  m_ngluons = NumberOfGluons();
  Flavour flav   = Flavour(kf_gluon);
  for (size_t i=0;i<m_ngluons;i++)  m_flavours.push_back(flav);
  for (size_t i=1;i<6;i++) {
    flav = Flavour(i);
    if (flav.Mass(true)>m_threshold) continue;
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
  for (size_t i=1;i<6;i++) {
    flav = Flavour(i);
    if (flav.Mass(true)>m_threshold) continue;
    totmass += 2.*flav.Mass(true);
    if (totmass>m_Ehat) break;
    m_nquarks+=2;
    if (flav.Bar()!=m_flavs[0] && flav.Bar()!=m_flavs[1]) m_flavours.push_back(flav);
    flav = flav.Bar();
    if (flav.Bar()!=m_flavs[0] && flav.Bar()!=m_flavs[1]) {
      m_flavours.push_back(flav);
      m_nquarks++;
      totmass += flav.Mass(true);
    }
  }
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
  size_t nquarkpairs = (m_flavours.size()-m_ngluons-2)/2., ncolours = nquarkpairs+m_ngluons+2;
  std::vector<size_t> cols[2];
  for (size_t j=0;j<2;j++) {
    for (size_t i=0;i<ncolours;i++) 
      cols[j].push_back(500+i);
  }
  for (size_t i=0;i<m_colours.size();i++) m_colours[i].clear(); m_colours.clear();
  m_colours.resize(m_flavours.size());
  size_t  pos, col, index;
  Flavour flav;
  for (size_t i=0;i<m_flavours.size()-1;++i) {
    flav = m_flavours[i];
    m_colours[i].resize(2);
    if ((flav.IsQuark() && !flav.IsAnti()) || flav.IsGluon()) {
      index = i<2?1:0;
      do {
	pos = size_t(cols[index].size()*ran->Get());
      } while (pos>=cols[index].size());
      m_colours[i][0] = cols[index][pos];
      if (cols[index].size()>1) {
	for (size_t j=pos;j<cols[index].size()-1;j++) {
	  cols[index][j] = cols[index][j+1];
	}
      }
      cols[index].pop_back();
    }
    else m_colours[i][0] = 0;
    if ((flav.IsQuark() && flav.IsAnti()) || flav.IsGluon()) {
      index = i<2?0:1;
      do {
	pos = size_t(cols[index].size()*ran->Get());
	if (pos>=cols[index].size()) continue;
        col = cols[index][pos];
      } while (col==m_colours[i][0]);
      m_colours[i][1] = col;
      if (cols[index].size()>1) {
	for (size_t j=pos;j<cols[index].size()-1;j++) {
	  cols[index][j] = cols[index][j+1];
	}
      }
      cols[index].pop_back();
    }
    else m_colours[i][1] = 0;
  }
  flav = m_flavours[m_flavours.size()-1];
  m_colours[m_flavours.size()-1].resize(2);
  m_colours[m_flavours.size()-1][0] = (((flav.IsQuark() && !flav.IsAnti()) || flav.IsGluon())?
				       cols[0].back():0);
  m_colours[m_flavours.size()-1][1] = (((flav.IsQuark() && flav.IsAnti()) || flav.IsGluon())?
				       cols[1].back():0);  
  if (p_colours[m_flavours.size()-1][0]==p_colours[m_flavours.size()-1][1]) {
    size_t help = p_colours[m_flavours.size()-1][1];
    p_colours[m_flavours.size()-1][1] = p_colours[2][1];
    p_colours[2][1] = help;
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
  exit(1);
}

bool XS_instanton::SetColours(const Vec4D_Vector& mom) {
  msg_Error()<<METHOD<<" shouldn't be called. Will terminate run.\n";
  exit(1);
}

DECLARE_TREEME2_GETTER(XS_instanton,"XS_instanton")
Tree_ME2_Base *ATOOLS::Getter<Tree_ME2_Base,External_ME_Args,XS_instanton>::
operator()(const External_ME_Args& args) const
{
  const Flavour_Vector fl=args.Flavours();
  if (fl.size()!=3) return NULL;
  if (!(fl[0].Strong() && fl[1].Strong() &&
	fl[2]==Flavour(kf_instanton))) return NULL;
  return new XS_instanton(args);
}
