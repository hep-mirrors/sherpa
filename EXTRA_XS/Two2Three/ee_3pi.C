#include "PHASIC++/Process/External_ME_Args.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Settings.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "MODEL/UFO/UFO_Model.H"
#include "MODEL/Main/Model_Base.H"

#include "EXTRA_XS/Main/ME2_Base.H"
#include "METOOLS/HadronCurrents/V_0_Isoscalar3Pi.H"
#include "METOOLS/HadronCurrents/Tools.H"
#include "METOOLS/Main/XYZFuncs.H"

using namespace EXTRAXS;
using namespace MODEL;
using namespace ATOOLS;
using namespace PHASIC;
using namespace std;

namespace EXTRAXS {

  ///////////////////////////////////////////////////////////////////////////
  //
  //  e+ e-  ->  gamma*  ->  pi+ pi- pi0
  //
  //  The hadronic side is anomalous, so the current is fixed up to one scalar
  //  form factor:
  //
  //      H^mu = eps^{mu al be ga} p_+al p_-be p_0ga  F(s; s_+-, s_+0, s_-0)
  //
  //  Both sides are Current_Base objects and the matrix element is the
  //  helicity sum of their contraction,
  //
  //      <|M|^2> = (4 pi alpha)^2 / s^2  1/4 Sum_{h1 h2} | L_{h1 h2} . H |^2 ,
  //
  //  which is what Current_ME does for the tau.  The leptonic current comes
  //  from the same XYZFunc spinors HADRONS++ uses, so the electron mass is
  //  carried exactly rather than dropped as it is in the massless tensor
  //  L^{mu nu} = 4 (k1^mu k2^nu + k2^mu k1^nu - g^{mu nu} k1.k2).  H is
  //  orthogonal to q = k1 + k2 by construction, so no gauge term survives.
  //
  //  NORMALISATION.  The gamma-3pi coupling is not perturbative, so the
  //  absolute scale of this process has to be fitted whatever we do.  It is
  //  collected into one factor, EE3Pi_norm, which also absorbs the spin-average
  //  convention and units.  What is derived from first principles here is the
  //  SHAPE -- the s dependence and the Dalitz structure -- and that is what the
  //  fit is then testing.
  //
  ///////////////////////////////////////////////////////////////////////////

  class ee_3pi : public ME2_Base {
  private:
    METOOLS::V_0_Isoscalar3Pi * p_had;   // gamma* -> 3 pi
    METOOLS::XYZFunc          * p_lept;  // e+ e- -> gamma*
    METOOLS::GeneralModel       m_model;
    Flavour_Vector              m_ffflavs;
    std::vector<int>            m_lepidx;
    double                       m_norm, m_alpha;
    int                          m_ip, m_im, m_i0;

    void ReadParameters();
  public:
    ee_3pi(const External_ME_Args & args);
    ~ee_3pi();
    double operator()(const ATOOLS::Vec4D_Vector & mom);
  };
}

void ee_3pi::ReadParameters() {
  // The form factor reads its couplings out of a GeneralModel, which is just a
  // name -> double map, so the run card is copied into one here.  Same tags the
  // form factor documents: EE3Pi_amp_<res> and EE3Pi_phase_<res>.
  Scoped_Settings s{ Settings::GetMainSettings()["EE_3PI"] };
  // Only one model here: hep-ph/0512180 Eq. (10) over registry rho line shapes.
  // The factorised variants that lived alongside it upstream are not ported --
  // this one supersedes them (chi2/ndf 1.80 against their best of 2.56).
  m_norm = s["Norm"].SetDefault(1.).Get<double>();

  // A--F, refitted against Belle 2024 with these line shapes.  NOT the paper's
  // published values, which assume its own inline propagators.
  const string tags[6] = { "A", "B", "C", "D", "E", "F" };
  const double def[6]  = { 18.20, -0.87, -0.5785, -1.2062, -0.72, -0.3947 };
  for (size_t i(0);i<6;++i)
    m_model["EE3Pi_"+tags[i]] = s[tags[i]].SetDefault(def[i]).Get<double>();

  // Isoscalar poles.  These are MODEL parameters, not particle properties --
  // see FF_0_Isoscalar3Pi.H for why.  Zero is not special here; the form
  // factor supplies the defaults.
  for (auto k : { "M_omega", "G_omega", "M_phi", "G_phi",
                  "M_omegaP", "G_omegaP", "M_omegaPP", "G_omegaPP" }) {
    Scoped_Settings v = s[k];
    const double d = v.SetDefault(-1.).Get<double>();
    if (d>0.) m_model[string("EE3Pi_")+k] = d;
  }

  // The first four are the isoscalar tower in s; the last two are the rho
  // tower in the two-pion sub-masses.  Defaults are documented where the form
  // factor reads them, in FF_0_Isoscalar3Pi::FixParameters.
  msg_Out()<<METHOD<<": hep-ph/0512180 Eq.(10) over registry rho line shapes"
	   <<", norm = "<<m_norm<<", A-F =";
  for (size_t i(0);i<6;++i) msg_Out()<<" "<<m_model["EE3Pi_"+tags[i]];
  msg_Out()<<"\n";
}

ee_3pi::ee_3pi(const External_ME_Args & args) :
  ME2_Base(args), p_had(NULL), p_lept(NULL), m_norm(1.), m_alpha(0.),
  m_ip(-1), m_im(-1), m_i0(-1)
{
  // 0, not 1: the S1/T1/U1 channels EXTRA_XS builds from this mask are all
  // 2 -> 2 only, and S1Channel throws outright for nout = 3.
  m_sintt = 0;
  m_oew   = 0;   // the coupling is carried explicitly below, not by the order
  m_oqcd  = 0;

  m_alpha = MODEL::s_model->ScalarConstant("alpha_QED");

  // m_flavs, not m_flavours: the latter is declared in ME2_Base but never
  // filled by its constructor.
  for (size_t i(2);i<m_flavs.size();++i) {
    const Flavour & fl = m_flavs[i];
    if      (fl.Kfcode()==kf_pi_plus && fl.IntCharge()>0) m_ip = i;
    else if (fl.Kfcode()==kf_pi_plus && fl.IntCharge()<0) m_im = i;
    else if (fl.Kfcode()==kf_pi)                          m_i0 = i;
  }
  if (m_ip<0 || m_im<0 || m_i0<0)
    THROW(fatal_error,"ee_3pi: final state is not pi+ pi- pi0.");

  ReadParameters();

  // Current_Base keeps a REFERENCE to the flavours, so they have to outlive
  // it -- hence the member copy rather than m_flavs, which ME2_Base owns.
  m_ffflavs = m_flavs;
  const vector<int> had = { m_ip, m_im, m_i0 };
  p_had = new METOOLS::V_0_Isoscalar3Pi(m_ffflavs,had,"V_0_Isoscalar3Pi");
  p_had->SetModelParameters(m_model);

  // Index 0 of the leptonic current is the barred spinor, i.e. the positron.
  // XYZFunc keeps a POINTER into this vector rather than copying it, so it
  // has to be a member and not a local.
  m_lepidx = { 1, 0 };
  p_lept = new METOOLS::XYZFunc(m_ffflavs,m_lepidx);
}

ee_3pi::~ee_3pi() {
  if (p_had)  { delete p_had;  p_had  = NULL; }
  if (p_lept) { delete p_lept; p_lept = NULL; }
}

double ee_3pi::operator()(const ATOOLS::Vec4D_Vector & mom) {
  const double s = (mom[0]+mom[1]).Abs2();
  if (s<=0.) return 0.;

  // The hadronic current, epsilon structure and all.  It is built by the
  // isobar recursion from declared spins and orbital angular momenta, not
  // written out here.
  p_had->Calc(mom,false);
  const Vec4C H = p_had->Get(size_t(0));

  // The leptonic current, one Vec4C per helicity pair.  A photon is pure
  // vector, so the right and left couplings are equal; their common size is
  // the explicit e^2 below.
  p_lept->Prepare(mom,false);
  const Complex one(1.,0.);
  double sum(0.);
  for (int h0(0);h0<2;++h0) {
    for (int h1(0);h1<2;++h1) {
      sum += std::norm(p_lept->L(0,h0,1,h1,one,one)*H);
    }
  }
  if (sum<=0.) return 0.;   // outside the physical region

  // 1/4 averages over the two incoming helicities each.
  return m_norm * sqr(4.*M_PI*m_alpha)/sqr(s) * 0.25*sum;
}

DECLARE_TREEME2_GETTER(EXTRAXS::ee_3pi,"ee_3pi")
Tree_ME2_Base *ATOOLS::Getter<PHASIC::Tree_ME2_Base,PHASIC::External_ME_Args,
			      EXTRAXS::ee_3pi>::
operator()(const External_ME_Args &args) const
{
  if (dynamic_cast<UFO::UFO_Model*>(MODEL::s_model)) return NULL;

  const Flavour_Vector fl = args.Flavours();
  if (fl.size()!=5) return NULL;
  if (!(fl[0]==Flavour(kf_e) && fl[1]==fl[0].Bar())) return NULL;

  // exactly one pi+, one pi- and one pi0, in any order
  int np(0), nm(0), nz(0);
  for (size_t i(2);i<5;++i) {
    if      (fl[i].Kfcode()==kf_pi_plus && fl[i].IntCharge()>0) ++np;
    else if (fl[i].Kfcode()==kf_pi_plus && fl[i].IntCharge()<0) ++nm;
    else if (fl[i].Kfcode()==kf_pi)                             ++nz;
    else return NULL;
  }
  if (np!=1 || nm!=1 || nz!=1) return NULL;

  return new ee_3pi(args);
}
