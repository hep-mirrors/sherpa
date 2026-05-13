#ifndef HADRON_RESCATTERING_XSecs_MesonMeson_H
#define HADRON_RESCATTERING_XSecs_MesonMeson_H

#include "HADRON_RESCATTERING/XSecs/HR_Parameters.H"
#include "HADRON_RESCATTERING/XSecs/ScatteringBase.H"
#include "HADRON_RESCATTERING/XSecs/Helper.H"
#include <map>
#include <list>
#include <tuple>
#include <vector>

namespace HADRON_RESCATTERING{

class MesonMeson : public ScatteringBase {
public:
    MesonMeson();
    ~MesonMeson();

    double XStot(const ATOOLS::Flavour& A, const ATOOLS::Flavour& B, const double& s) override;
    double XSel (const ATOOLS::Flavour& A, const ATOOLS::Flavour& B, const double& s);

    // Component accessors (kept for Tests/plotting)
    std::vector<double> PiPlusPiMinusResComponents(double const& s);
    std::vector<double> PiPlusPiZeroResComponents (double const& s);
    std::vector<double> PiZeroPiZeroResComponents (double const& s);
    std::vector<double> KPlusPiMinusResComponents (double const& s);

private:
    struct ChannelConfig {
        const std::vector<ResonanceParams>* resonances;
        double m1;
        double m2;
        double threshold;
        bool   addF500;   // true for pi+pi- and pi0pi0
    };

    using ChannelKey = std::tuple<int, int, bool, bool>; 
    // (kfA, kfB, AisAnti, BisAnti) — both are mesons so no baryon/meson split
    std::map<ChannelKey, ChannelConfig> m_channels;

    void buildChannelMap();

    double computeChannel(const double& s, const ChannelConfig& cfg);

    double calculateResonanceSigma(double sqrt_s, const ResonanceParams& res,
                                   double m1, double m2) const;

    std::vector<double> resonanceComponents(const double& s,
                                            const ChannelConfig& cfg) const;

    // f500 special low-energy table — kept separate
    double PiPlusPiMinusf500(double const& s);

    bool m_test;
    void Tests();
};

}

#endif

