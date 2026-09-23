#include "METOOLS/HadronCurrents/FormFactors/FF_PP_Wave.H"
#include "METOOLS/HadronCurrents/FormFactors/Line_Shapes.H"
#include "METOOLS/HadronCurrents/FormFactors/Resonance_Base.H"
#include "METOOLS/HadronCurrents/Tools.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"

using namespace METOOLS;
using namespace ATOOLS;
using namespace std;

FF_PP_Wave::FF_PP_Wave(const FF_Parameters & params) :
  FormFactor_Base(params),
  m_model(pp_wave_model::LASS), m_L(0),
  m_m1(0.), m_m2(0.), m_sthr(0.),
  p_prop(NULL), p_lass(NULL),
  m_a(1.94), m_r(1.76), m_cutoff(-1.)
{
  m_m1 = m_masses[m_pi[0]];
  m_m2 = m_masses[m_pi[1]];
  m_sthr = sqr(m_m1+m_m2);
  ReadParameters();
}

FF_PP_Wave::~FF_PP_Wave() {
  if (p_prop) delete p_prop;
  if (p_lass) delete p_lass;
}

//////////////////////////////////////////////////////////////////////////////
// Propagator assembly
//////////////////////////////////////////////////////////////////////////////

Propagator_Base * FF_PP_Wave::BuildOne(const std::string & tag) {
  const GeneralModel & md = *p_model;
  const int  kf   = int(md(m_name+"_KF"+tag,0.));
  const int  type = int(md(m_name+"_TYPE"+tag,
                           double(int(resonance_type::running))));

  if (resonance_type(type)==resonance_type::flatte) {
    // Two-channel Flatte. Channel 1 is the one this current produces,
    // channel 2 the one that opens nearby - K Kbar for both f_0(980)
    // and a_0(980). Channel 1 defaults to the two pseudoscalars of
    // THIS final state, which is almost always what is wanted.
    const double M  = md(m_name+"_MFLAT"+tag,0.);
    const double g1 = md(m_name+"_G1"+tag,0.);
    const double g2 = md(m_name+"_G2"+tag,0.);
    if (M<=0. || g1<=0.)
      THROW(fatal_error,"Flatte requested for "+m_name+" but "+m_name+
            "_MFLAT"+tag+" / "+m_name+"_G1"+tag+" are not set. A Flatte "
            "mass and its couplings are correlated and only meaningful "
            "as a set - do not mix them with PDG Breit-Wigner values.");
    const double ma1 = md(m_name+"_MA1"+tag,m_m1);
    const double mb1 = md(m_name+"_MB1"+tag,m_m2);
    const double ma2 = md(m_name+"_MA2"+tag,0.493677);   // K+
    const double mb2 = md(m_name+"_MB2"+tag,0.493677);
    Flavour fl(kf_none);
    if (kf!=0) { fl = Flavour((kf_code)abs(kf)); if (kf<0) fl = fl.Bar(); }
    return new Flatte(M,g1,g2,ma1,mb1,ma2,mb2,fl);
  }

  if (resonance_type(type)==resonance_type::complex_pole) {
    // No Total_Width_Base is consulted: a complex pole has no running
    // width by construction. The pole position must be given.
    const double Mp = md(m_name+"_MPOLE"+tag,0.);
    const double Gp = md(m_name+"_GPOLE"+tag,0.);
    if (Mp<=0. || Gp<=0.)
      THROW(fatal_error,"Complex pole "+m_name+"_TYPE"+tag+" requested but "
            +m_name+"_MPOLE"+tag+" / "+m_name+"_GPOLE"+tag+" are not set. "
            "A complex pole cannot fall back on a Breit-Wigner mass: it is "
            "used precisely where no real-axis mass exists.");
    // Carry the kf if one was given, purely so the diagnostic dump
    // names the state instead of printing "no_particle".
    Flavour fl(kf_none);
    if (kf!=0) { fl = Flavour((kf_code)abs(kf)); if (kf<0) fl = fl.Bar(); }
    return new Complex_Pole(Mp,Gp,fl);
  }

  if (kf==0)
    THROW(fatal_error,"No resonance code given for "+m_name+"_KF"+tag+".");
  Flavour fl((kf_code)abs(kf));
  if (kf<0) fl = fl.Bar();
  Total_Width_Base * w = LineShapes->Get(fl);
  if (w==NULL)
    THROW(fatal_error,"No lineshape registered for "+ToString(fl)+
          " (requested by "+m_name+"_KF"+tag+"). Add it in "
          "Line_Shapes::Init() rather than parametrising it here - the "
          "running width must be the same object the rest of the library "
          "uses for this resonance.");
  // Gounaris-Sakurai needs the dominant daughter mass; for the pi pi P
  // wave that is the pion, which is m_m1 here.
  if (resonance_type(type)==resonance_type::GS)
    return new BreitWigner(w,resonance_type::GS,m_m1);
  return new BreitWigner(w,resonance_type(type));
}

void FF_PP_Wave::ReadParameters() {
  if (p_model==NULL) return;
  const GeneralModel & md = *p_model;
  const string & n = m_name;

  m_model = pp_wave_model(int(md(n+"_MODEL",
                                 double(int(pp_wave_model::LASS)))));
  m_L     = int(md(n+"_L",0.));
  if (m_L<0 || m_L>2)
    THROW(fatal_error,"FF_PP_Wave '"+n+"': only L = 0, 1, 2 implemented.");
  if (m_L!=0 && m_model==pp_wave_model::LASS)
    THROW(fatal_error,"FF_PP_Wave '"+n+"': LASS is an S-wave "
          "parametrization but L = "+ToString(m_L)+" was requested.");

  m_a      = md(n+"_LASS_A",m_a);
  m_r      = md(n+"_LASS_R",m_r);
  m_cutoff = md(n+"_CUTOFF",-1.);

  switch (m_model) {
  case pp_wave_model::flat:
    break;
  case pp_wave_model::resonance:
    p_prop = BuildOne("0");
    break;
  case pp_wave_model::two_resonance: {
    // Assembled with the library's own Summed_Propagator, so that the
    // admixture behaves exactly like the rho/rho'/rho'' and
    // kappa/K*_0(1430) towers already built in FF_0_PP.
    Summed_Propagator * sum = new Summed_Propagator();
    sum->Add(BuildOne("0"),Complex(1.,0.));
    const double amp = md(n+"_AMP1",1.), ph = md(n+"_PHASE1",0.);
    sum->Add(BuildOne("1"),Complex(amp*cos(ph),amp*sin(ph)));
    p_prop = sum;
    break;
  }
  case pp_wave_model::LASS:
    p_lass = BuildOne("0");
    break;
  default:
    THROW(fatal_error,"FF_PP_Wave '"+n+"': unknown model "+ToString(int(m_model)));
  }

  DumpPropagatorStructure(n,int(m_model),
                          (p_prop!=NULL ? p_prop : p_lass));
  PrintSetup();
  // p_model refers to a by-value GeneralModel in the owning current's
  // SetModelParameters frame and dangles from here on - see FF_P_X.
  p_model = NULL;
}

//////////////////////////////////////////////////////////////////////////////

double FF_PP_Wave::Breakup(const double & s) const {
  if (s<=m_sthr) return 0.;
  const double num = (s-sqr(m_m1+m_m2))*(s-sqr(m_m1-m_m2));
  return (num>0. ? sqrt(num)/(2.*sqrt(s)) : 0.);
}

Complex FF_PP_Wave::LASSAmp(const double & s) {
  // A_LASS = sin(dB) e^{i dB} + e^{2 i dB} D(s), with the elastic
  // effective-range background cot dB = 1/(a p) + r p / 2. The first
  // term is the non-resonant piece that dominates just above threshold
  // and that a bare resonance pole cannot reproduce.
  const double p = Breakup(s);
  if (p<1.e-12) return Complex(0.,0.);
  const double cotdB = 1./(m_a*p) + 0.5*m_r*p;
  const double dB    = atan2(1.,cotdB);
  const Complex ei(cos(dB),sin(dB));
  return sin(dB)*ei + ei*ei*(*p_lass)(s);
}

Complex FF_PP_Wave::Value(const double & s) {
  double se = s;
  if (m_cutoff>0. && s>m_cutoff) se = m_cutoff;
  if (se<=m_sthr) return Complex(0.,0.);
  switch (m_model) {
  case pp_wave_model::flat:          return Complex(1.,0.);
  case pp_wave_model::resonance:
  case pp_wave_model::two_resonance: return (*p_prop)(se);
  case pp_wave_model::LASS:          return LASSAmp(se);
  default: break;
  }
  return Complex(0.,0.);
}

Complex FF_PP_Wave::operator()(const ATOOLS::Vec4D_Vector& moms) {
  // s from the SUM of the two pseudoscalars - the FF_0_PP convention.
  return Value((moms[m_pi[0]]+moms[m_pi[1]]).Abs2());
}

void FF_PP_Wave::PrintSetup() const {
  msg_Tracking()<<"### FF_PP_Wave '"<<m_name<<"': L = "<<m_L
                <<", model "<<int(m_model)
                <<", threshold sqrt(s) = "<<sqrt(m_sthr)<<" GeV\n";
  if (m_model==pp_wave_model::LASS)
    msg_Tracking()<<"###   LASS background a = "<<m_a<<" GeV^-1, r = "
                  <<m_r<<" GeV^-1\n";
  msg_Tracking()<<"###   running widths taken from the LineShapes "
                <<"registry; complex-pole and Flatte propagators carry "
                <<"their own parameters.\n";
}

DECLARE_FF_GETTER(FF_PP_Wave,"FF_PP_Wave")

FormFactor_Base * ATOOLS::Getter<FormFactor_Base,FF_Parameters,FF_PP_Wave>::
operator()(const METOOLS::FF_Parameters &params) const
{
  if (params.m_pi.size()!=2) return NULL;
  return new FF_PP_Wave(params);
}
