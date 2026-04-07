#include "HADRON_RESCATTERING/XSecs/MesonMeson.H"
#include "ATOOLS/Math/Histogram.H"
#include "ATOOLS/Org/Message.H"
#include <algorithm>

using namespace HADRON_RESCATTERING;
using namespace ATOOLS;
using namespace std;


MesonMeson::MesonMeson() : ScatteringBase()
{
    m_test = false;
    MesonMeson::buildChannelMap();
    if (m_test) {
        Tests();
        msg_Out() << "Initializing MesonMeson ctor" << std::endl;
    }
}

MesonMeson::~MesonMeson() {} 


void MesonMeson::buildChannelMap()
{
    //ADD ALL ORDERINGS.

    // pi+ pi-  (one anti, one not — both orderings)
    m_channels[{211, 211, false, true}]  = { &m_PiPlusPiMinusRes, PION_MASS, PION_MASS, 2.*PION_MASS, true  };
    m_channels[{211, 211, true, false}]  = { &m_PiPlusPiMinusRes, PION_MASS, PION_MASS, 2.*PION_MASS, true  };

    // pi+ pi+  and pi- pi- ADD THEM

    // pi+/- pi0  (all four orderings)
    m_channels[{211, 111, false, false}] = { &m_PiPlusPiZeroRes, PION_MASS, PION_MASS, 2.*PION_MASS, false };
    m_channels[{211, 111, true,  false}] = { &m_PiPlusPiZeroRes, PION_MASS, PION_MASS, 2.*PION_MASS, false };
    m_channels[{111, 211, false, false}] = { &m_PiPlusPiZeroRes, PION_MASS, PION_MASS, 2.*PION_MASS, false };
    m_channels[{111, 211, false, true }] = { &m_PiPlusPiZeroRes, PION_MASS, PION_MASS, 2.*PION_MASS, false };

    // pi0 pi0
    m_channels[{111, 111, false, false}] = { &m_PiZeroPiZeroRes, PION_MASS, PION_MASS, 2.*PION_MASS, true  };

    // K+ pi-  and pi- K+  (both orderings)
    m_channels[{321, 211, false, true }] = { &m_KPlusPiMinusRes, CHARGEDKAON_MASS, PION_MASS, CHARGEDKAON_MASS+PION_MASS, false };
    m_channels[{211, 321, true,  false}] = { &m_KPlusPiMinusRes, PION_MASS, CHARGEDKAON_MASS, CHARGEDKAON_MASS+PION_MASS, false };
}

double MesonMeson::XStot(const Flavour& A, const Flavour& B, const double& s)
{
    if (!(A.IsMeson() && B.IsMeson())) return 0.;

    ChannelKey key { A.Kfcode(), B.Kfcode(), A.IsAnti(), B.IsAnti() };

    std::map<ChannelKey, ChannelConfig>::iterator it = m_channels.find(key);
    if (it != m_channels.end())
        return computeChannel(s, it->second);

    msg_Out() << "MesonMeson::XStot: no channel for ("
              << A.Kfcode() << (A.IsAnti()?"~":"") << ", "
              << B.Kfcode() << (B.IsAnti()?"~":"") << ")"
              << " isAnti: " << A.IsAnti() << " " << B.IsAnti() << std::endl;
    return 0.;
}

double MesonMeson::XSel(const Flavour&, const Flavour&, const double&)
{
    return 0.; //ADD THIS FUNCTION
}


double MesonMeson::computeChannel(const double& s, const ChannelConfig& cfg)
{
    if (s <= 0.) return 0.;
    double sqrt_s = std::sqrt(s);
    if (sqrt_s <= cfg.threshold) return 0.;

    double total = 0.;
    for (const ResonanceParams& res : *cfg.resonances)
        total += calculateResonanceSigma(sqrt_s, res, cfg.m1, cfg.m2);

    double result = total * GevToMB;
    if (cfg.addF500) result += PiPlusPiMinusf500(s);
    return result;
}


double MesonMeson::calculateResonanceSigma(double sqrt_s, const ResonanceParams& res, double m1, double m2) const
{
    double pcm = p_CM(sqrt_s, m1, m2);
    if (pcm <= 0.) return 0.;

    double Gamma_tot = totalWidth(sqrt_s, res.mass, res.width0, res.branchingRatio,
                                  res.L, res.m1, res.m2, res.LFactors);
    if (Gamma_tot <= 0.) return 0.;

    double spinFactor = res.twoSRplusOne;
    double massDiff   = res.mass - sqrt_s;
    double denom      = massDiff*massDiff + 0.25*Gamma_tot*Gamma_tot;
    if (denom <= 0.) return 0.;

    return spinFactor
         * (M_PI / (pcm*pcm))
         * (Gamma_tot * Gamma_tot) / denom
         * res.probFactor;
}

//COMPONENT VECTORS FOR PLOTTING PURPOSES ONLY, NOT USED IN XStot CALCULATION. NOT SURE IF IT WORKS CORRECTLY.
//come to think of it, it should not work because when called, I am only passing one ordering. 
//like PiPlusPiMinusResComponents only gets the ordering with pi+ as A and pi- as B, but not the other way around.
//Again, come to think of it, it should be fine because the resonances are the same for both orderings, and the calculateResonanceSigma function does not depend on the ordering. So it should give the same result for both orderings, which is what we want for plotting purposes.
std::vector<double> MesonMeson::resonanceComponents(const double& s, const ChannelConfig& cfg) const
{
    double sqrt_s = std::sqrt(s);
    std::vector<double> out;
    out.reserve(cfg.resonances->size());

    if (sqrt_s <= cfg.threshold) {
        out.assign(cfg.resonances->size(), 0.);
        return out;
    }

    for (const ResonanceParams& res : *cfg.resonances)
        out.push_back(calculateResonanceSigma(sqrt_s, res, cfg.m1, cfg.m2) * GevToMB);

    return out;
}

std::vector<double> MesonMeson::PiPlusPiMinusResComponents(double const& s)
{
    return resonanceComponents(s, m_channels.at({211, 211, false, true}));
}

std::vector<double> MesonMeson::PiPlusPiZeroResComponents(double const& s)
{
    return resonanceComponents(s, m_channels.at({211, 111, false, false}));
}

std::vector<double> MesonMeson::PiZeroPiZeroResComponents(double const& s)
{
    return resonanceComponents(s, m_channels.at({111, 111, false, false}));
}

std::vector<double> MesonMeson::KPlusPiMinusResComponents(double const& s)
{
    return resonanceComponents(s, m_channels.at({321, 211, false, true}));
}


double MesonMeson::PiPlusPiMinusf500(double const& s)
{
    double sqs = std::sqrt(s);

    static const std::vector<double> T = {
        0.279150, 0.293861, 0.308572, 0.323283, 0.337994,
        0.352706, 0.367417, 0.382128, 0.396839, 0.411550,
        0.426261, 0.440972, 0.455683, 0.470394, 0.485106,
        0.499817, 0.514528, 0.529239, 0.543950, 0.558661,
        0.573372, 0.588083, 0.602794, 0.617506, 0.632217,
        0.646928, 0.661639, 0.676350, 0.691061, 0.705772,
        0.720483, 0.735194, 0.749906, 0.764617, 0.779328,
        0.794039, 0.808750, 0.823461, 0.838172, 0.852883,
        0.867594, 0.882306, 0.897017, 0.911728, 0.926439,
        0.941150, 0.955861, 0.970572, 0.985283, 1.000000
    };

    static const std::vector<double> sigma_mbf = {
        8.15994,  9.53565, 11.0102,  12.5738,  14.2131,
        15.9117,  17.6494, 19.4033,  21.1478,  22.8556,
        24.4988,  26.0502, 27.4844,  28.7791,  29.9161,
        30.8821,  31.669,  32.2741,  32.6997,  32.9524,
        33.0426,  32.9833, 32.7894,  32.4769,  32.0621,
        31.5607,  30.9881, 30.3583,  29.6841,  28.9768,
        28.2464,  27.5015, 26.749,   25.995,   25.244,
        24.4994,  23.7637, 23.0379,  22.3219,  21.6142,
        20.9102,  20.182,  19.3792,  18.4299,  17.2295,
        15.6266,  13.4101, 10.3267,   6.29745,  3.39533
    };

    if (sqs <= T.front() || sqs >= T.back()) return 0.;

    std::vector<double>::const_iterator it = 
        std::upper_bound(T.begin(), T.end(), sqs);
    size_t i = std::distance(T.begin(), it) - 1;
    double t = (sqs - T[i]) / (T[i+1] - T[i]);
    return sigma_mbf[i] + t*(sigma_mbf[i+1] - sigma_mbf[i]);
}


void MesonMeson::Tests()
{
    msg_Out() << "MesonMeson::Tests called" << std::endl;

    size_t bins = 10000;
    double pmin = 0., pmax = 60., pinc = (pmax-pmin)/double(bins);

    map<string, Histogram*> histos;

    histos["PiPlusPiMinus_total"]    = new Histogram(0, pmin, pmax, bins);
    histos["PiPlusPi0_total"]        = new Histogram(0, pmin, pmax, bins);
    histos["Pi0Pi0_total"]           = new Histogram(0, pmin, pmax, bins);
    histos["KPlusPiMinus_total"]     = new Histogram(0, pmin, pmax, bins);

    histos["PiPlusPiMinusomegares"]   = new Histogram(0, pmin, pmax, bins);
    histos["PiPlusPiMinusrho770res"]  = new Histogram(0, pmin, pmax, bins);
    histos["PiPlusPiMinusf0980res"]   = new Histogram(0, pmin, pmax, bins);
    histos["PiPlusPiMinusf21270res"]  = new Histogram(0, pmin, pmax, bins);
    histos["PiPlusPiMinusf01370res"]  = new Histogram(0, pmin, pmax, bins);
    histos["PiPlusPiMinusf21525res"]  = new Histogram(0, pmin, pmax, bins);
    histos["PiPlusPiMinusrho1465res"] = new Histogram(0, pmin, pmax, bins);
    histos["PiPlusPiMinusrho1700res"] = new Histogram(0, pmin, pmax, bins);

    histos["PiPlusPiZerorho770res"]  = new Histogram(0, pmin, pmax, bins);
    histos["PiPlusPiZerorho1465res"] = new Histogram(0, pmin, pmax, bins);
    histos["PiPlusPiZerorho1700res"] = new Histogram(0, pmin, pmax, bins);

    histos["PiZeroPiZeroomegares"]  = new Histogram(0, pmin, pmax, bins);
    histos["PiZeroPiZerof0980res"]  = new Histogram(0, pmin, pmax, bins);
    histos["PiZeroPiZerof01370res"] = new Histogram(0, pmin, pmax, bins);
    histos["PiZeroPiZerof21270res"] = new Histogram(0, pmin, pmax, bins);
    histos["PiZeroPiZerof21525res"] = new Histogram(0, pmin, pmax, bins);

    histos["KPlusPiMinusK890res"]     = new Histogram(0, pmin, pmax, bins);
    histos["KPlusPiMinusKstar0res"]   = new Histogram(0, pmin, pmax, bins);
    histos["KPlusPiMinusKstar1430res"]= new Histogram(0, pmin, pmax, bins);
    histos["KPlusPiMinusKstar1410res"]= new Histogram(0, pmin, pmax, bins);
    histos["KPlusPiMinusKstar1680res"]= new Histogram(0, pmin, pmax, bins);

    Flavour flA(kf_pi_plus); flA = flA.Bar();
    Flavour flB(kf_K_plus);

    for (int i = 0; i < (int)bins; ++i) {
        double plab = pmin + i*pinc;
        double s    = sqr(flA.HadMass()) + sqr(flB.HadMass())
                    + 2.*flA.HadMass()*sqrt(sqr(flB.HadMass()) + sqr(plab));

        histos["PiPlusPiMinus_total"]->Insert(i, computeChannel(s, m_channels.at({211,211,false,true})));
        histos["PiPlusPi0_total"]    ->Insert(i, computeChannel(s, m_channels.at({211,111,false,false})));
        histos["Pi0Pi0_total"]       ->Insert(i, computeChannel(s, m_channels.at({111,111,false,false})));
        histos["KPlusPiMinus_total"] ->Insert(i, computeChannel(s, m_channels.at({321,211,false,true})));

        std::vector<double> pm = PiPlusPiMinusResComponents(s);
        histos["PiPlusPiMinusomegares"]  ->Insert(i, pm[0]);
        histos["PiPlusPiMinusrho770res"] ->Insert(i, pm[1]);
        histos["PiPlusPiMinusf0980res"]  ->Insert(i, pm[2]);
        histos["PiPlusPiMinusf21270res"] ->Insert(i, pm[3]);
        histos["PiPlusPiMinusf01370res"] ->Insert(i, pm[4]);
        histos["PiPlusPiMinusf21525res"] ->Insert(i, pm[5]);
        histos["PiPlusPiMinusrho1465res"]->Insert(i, pm[6]);
        histos["PiPlusPiMinusrho1700res"]->Insert(i, pm[7]);

        std::vector<double> pz = PiPlusPiZeroResComponents(s);
        histos["PiPlusPiZerorho770res"] ->Insert(i, pz[0]);
        histos["PiPlusPiZerorho1465res"]->Insert(i, pz[1]);
        histos["PiPlusPiZerorho1700res"]->Insert(i, pz[2]);

        std::vector<double> zz = PiZeroPiZeroResComponents(s);
        histos["PiZeroPiZeroomegares"] ->Insert(i, zz[0]);
        histos["PiZeroPiZerof0980res"] ->Insert(i, zz[1]);
        histos["PiZeroPiZerof01370res"]->Insert(i, zz[2]);
        histos["PiZeroPiZerof21270res"]->Insert(i, zz[3]);
        histos["PiZeroPiZerof21525res"]->Insert(i, zz[4]);

        std::vector<double> kp = KPlusPiMinusResComponents(s);
        histos["KPlusPiMinusK890res"]    ->Insert(i, kp[0]);
        histos["KPlusPiMinusKstar0res"]  ->Insert(i, kp[1]);
        histos["KPlusPiMinusKstar1430res"]->Insert(i, kp[2]);
        histos["KPlusPiMinusKstar1410res"]->Insert(i, kp[3]);
        histos["KPlusPiMinusKstar1680res"]->Insert(i, kp[4]);
    }

    for (std::map<string,Histogram*>::iterator hit = histos.begin();
         hit != histos.end(); ++hit) {
        string name = "XSecs_MesonMesonTest/" + hit->first + ".dat";
        hit->second->Output(name);
        delete hit->second;
    }
    histos.clear();
}