#include "HADRON_RESCATTERING/XSecs/BaryonMeson.H"
#include "HADRON_RESCATTERING/XSecs/HPR1R2.H"
#include "HADRON_RESCATTERING/XSecs/Helper.H"
#include "HADRON_RESCATTERING/XSecs/ScatteringBase.H"
#include "ATOOLS/Math/Histogram.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/MathTools.H"

using namespace HADRON_RESCATTERING;
using namespace ATOOLS;
using namespace std;

BaryonMeson::BaryonMeson()  : ScatteringBase()
{
    m_test = false;
    if (m_test) 
    {
        msg_Out() << "Baryon Meson ctor called."<< std::endl;
        BaryonMeson::Tests(); 
        exit(1); 
    }
    BaryonMeson::buildChannelMap();  // vectors already filled by ScatteringBase ctor
}

void BaryonMeson::buildChannelMap()
{
    // { baryon kfcode, meson kfcode, meson isAnti } -> ChannelConfig
    m_channels = {

        // ── Proton channels ──────────────────────────────────────────────────
        { {2212, 211, false}, { &pPiPlus_nPiMinus, PROTON_MASS,  PION_MASS,
                                hpr1r2::pPiPlus,  -1,
                                PROTON_MASS + PION_MASS, 2. } },

        { {2212, 211, true},  { &pPiMZ_nPiPZ,     PROTON_MASS,  PION_MASS,
                                hpr1r2::pPiMinus,   0,
                                PROTON_MASS + PION_MASS, 2. } },

        { {2212, 111, false}, { &pPiMZ_nPiPZ,     PROTON_MASS,  PION_MASS,
                                hpr1r2::pPiZero,    1,
                                PROTON_MASS + PION_MASS, 2. } },

        // ── Neutron channels ─────────────────────────────────────────────────
        { {2112, 211, false}, { &pPiMZ_nPiPZ,     NEUTRON_MASS, PION_MASS,
                                hpr1r2::nPiPlus,    0,
                                NEUTRON_MASS + PION_MASS, 2. } },

        { {2112, 211, true},  { &pPiPlus_nPiMinus, NEUTRON_MASS, PION_MASS,
                                hpr1r2::nPiMinus,  -1,
                                NEUTRON_MASS + PION_MASS, 2. } },

        { {2112, 111, false}, { &pPiMZ_nPiPZ,     NEUTRON_MASS, PION_MASS,
                                hpr1r2::nPiZero,    1,
                                NEUTRON_MASS + PION_MASS, 2. } },

        { {2112, 321, true},  { &m_Resonances_n_KaonMinus, NEUTRON_MASS, KAON_MASS,
                                hpr1r2::nKMinus,   -1,
                                NEUTRON_MASS + KAON_MASS, 2.16 } },
    };
}

BaryonMeson::~BaryonMeson() {}


double BaryonMeson::XStot(const Flavour& A, const Flavour& B, const double& s)
{
    if (!((A.IsMeson() && B.IsBaryon()) || (A.IsBaryon() && B.IsMeson())))
        return 0.;
    Flavour baryon ;// = A.IsBaryon() ? A : B;
    Flavour meson  ;//= B.IsMeson()  ? B : A;

    if (A.IsBaryon()) {baryon = A; meson  = B;}
    else              {baryon = B; meson  = A;}

    ChannelKey key { baryon.Kfcode(), meson.Kfcode(), meson.IsAnti() };

    std::map<ChannelKey, ChannelConfig>::iterator it = m_channels.find(key);
    if (it != m_channels.end())
        return computeChannel(s, it->second);

    // pK- kept separate: has its own low-energy resonance formula
    if (baryon.Kfcode() == 2212 && meson.Kfcode() == 321 && meson.IsAnti())
        return pKaonMinus(s);

    //msg_Out() << "BaryonMeson::XStot: no channel found for ("
    //        << baryon.Kfcode() << ", " << meson.Kfcode()
    //        << ", anti= " << meson.IsAnti() << ")" << std::endl;
    return 0.;
}

double BaryonMeson::computeChannel(const double& s, const ChannelConfig& cfg)
{
    if (s <= 0.)                       return 0.;
    double sqrt_s = std::sqrt(s);
    if (sqrt_s <= cfg.threshold)       return 0.;
    if (sqrt_s >= cfg.hpr1r2Threshold) return m_hpr1r2.xs_tot(cfg.hpr1r2Code, s);

    double total = 0.;
    for (const ResonanceParams& res  : *cfg.resonances)
        total += calculateResonanceSigma(sqrt_s, res, cfg.m1, cfg.m2, cfg.channelIndex);

    return total * GevToMB;
}

double BaryonMeson::calculateResonanceSigma(double sqrt_s, const ResonanceParams& res, double m1, double m2, int channelIndex) const
{
    double pcm = p_CM(sqrt_s, m1, m2);
    if (pcm <= 1e-3 || sqrt_s <= m1 + m2 + 1e-3) return 0.;

    double Gamma_tot = totalWidth(sqrt_s, res.mass, res.width0, res.branchingRatio,
                                  res.L, res.m1, res.m2, res.LFactors);
    if (Gamma_tot < 1e-9) return 0.;

    // Isospin weight depends on which sub-channel we are in
    double isospinWeight = 1.0;
    if      (channelIndex == 0) isospinWeight = 1. - res.probFactor;
    else if (channelIndex == 1) isospinWeight =      res.probFactor;

    double spinFactor = res.twoSRplusOne / 2.;
    double massDiff   = res.mass - sqrt_s;
    double denom      = massDiff*massDiff + 0.25*Gamma_tot*Gamma_tot;

    double sigma = isospinWeight
         * spinFactor
         * M_PI / (pcm * pcm)
         * (res.branchingRatio[0] * Gamma_tot * Gamma_tot) / denom;
    return sigma;
}


double BaryonMeson::pKaonMinus(const double& s)
{
    double sqrt_s = std::sqrt(s);
    if (sqrt_s <= PROTON_MASS + KAON_MASS) return 0.;

    if (sqrt_s >= 2.16)
        return m_hpr1r2.xs_tot(hpr1r2::pKMinus, s);

    double total = 0.;
    for (const ResonanceParams& res : m_Resonances)
        total += NKaonScatteringXSec(sqrt_s, res, PROTON_MASS, KAON_MASS);

    double e0 = 1.433;
    if      (sqrt_s < 1.4738188)
        total += 5.93763355 / pow(sqrt_s - 1.251377, 2);
    else if (sqrt_s < 1.485215)
        total += -1.296457765e7 * pow(sqrt_s - e0, 4)
               +  2.160975431e4 * pow(sqrt_s - e0, 2) + 120.;
    else if (sqrt_s < 1.977)
        total += 3. + 1.0777e6 * exp(-6.4463  * sqrt_s)
               - 10. * exp(-pow(sqrt_s - 1.644, 2) / 0.004)
               + 10. * exp(-pow(sqrt_s - 1.977, 2) / 0.004);
    else
        total += 1.0777e6 * exp(-6.44463 * sqrt_s) + 12.5;

    return total;   // already in mb
}

double BaryonMeson::NKaonScatteringXSec(double sqrt_s, const ResonanceParams& res, double m1, double m2) const
{
    double pcm = p_CM(sqrt_s, m1, m2);
    if (pcm <= 1e-3 || sqrt_s <= m1 + m2 + 1e-3) return 0.;

    double Gamma_tot = totalWidth(sqrt_s, res.mass, res.width0, res.branchingRatio,
                                  res.L, res.m1, res.m2, res.LFactors);
    if (Gamma_tot < 1e-9) return 0.;

    double spinFactor = res.twoSRplusOne / 2.;
    double massDiff   = res.mass - sqrt_s;
    double denom      = massDiff*massDiff + 0.25*Gamma_tot*Gamma_tot;

    return GevToMB * res.probFactor
         * spinFactor
         * M_PI / (pcm * pcm)
         * (res.branchingRatio[0] * Gamma_tot * Gamma_tot) / denom;
}

void BaryonMeson::Tests() {

  msg_Out() << "BaryonMeson::Test  called" << std::endl;
  size_t bins = 10000;
  double pmin = 0. , pmax = 30., pinc = (pmax-pmin)/double(bins);
  map<string,Histogram *>  histos;
  histos["protonPionPlus_Total"] = new Histogram(0,pmin,pmax,bins);
  Flavour flA(kf_p_plus);
  double  plab, s;

  Flavour flB(kf_pi_plus);

  for (int i=0;i<bins;i++) 
  {
    plab   = pmin+i*pinc;
    s      = ( sqr(flA.HadMass()) + sqr(flB.HadMass()) +
	       2.*flA.HadMass()*sqrt(sqr(flB.HadMass())+sqr(plab)) );
    histos["protonPionPlus_Total"]->Insert(i,   BaryonMeson::XStot(flA,flB,s));
  }

  Histogram * histo;
  std::string name;
  for (std::map<std::string,Histogram *>::iterator
	 hit=histos.begin();hit!=histos.end();hit++) {
    histo = hit->second;
    name  = std::string("XSecs_BaryonMesonTest/")+hit->first+std::string(".dat");
    histo->Output(name);
    delete histo;
  }
  histos.clear();
}



















