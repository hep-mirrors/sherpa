#ifndef HADRON_RESCATTERING_XSecs_BaryonMeson_H
#define HADRON_RESCATTERING_XSecs_BaryonMeson_H

#include "HADRON_RESCATTERING/XSecs/NNfits.H"
#include "HADRON_RESCATTERING/XSecs/HPR1R2.H"
#include "HADRON_RESCATTERING/XSecs/ScatteringBase.H"
#include "ATOOLS/Phys/Flavour.H"

namespace HADRON_RESCATTERING {
  class BaryonMeson  : public ScatteringBase 
  {
    private: 
    NNfits  m_BM;
    HPR1R2  m_hpr1r2;
    double  m_s, m_mA, m_mA2, m_mB, m_mB2;
    
    double xstot(long int & A, long int & B, const double & plab);
    double xsel(long int & A, long int & B, const double & plab);
    
    bool m_test = false;
    void Tests();

    struct ChannelConfig {
        const std::vector<ResonanceParams>* resonances;
        double       m1;
        double       m2;
        hpr1r2::code hpr1r2Code;
        int          channelIndex;     // -1 = plain BW, 0 = (1-probFactor), 1 = probFactor
        double       threshold;        // minimum sqrt_s (= m1 + m2)
        double       hpr1r2Threshold;  // sqrt_s above which we use HPR1R2
    };

    using ChannelKey = std::tuple<int, int, bool>; // (baryon kf, meson kf, meson isAnti)

    std::map<ChannelKey, ChannelConfig> m_channels;

    void   buildChannelMap();
    double computeChannel(const double& s, const ChannelConfig& cfg); 

  public:

    BaryonMeson();
    ~BaryonMeson();
  double calculateResonanceSigma(double sqrt_s, const ResonanceParams& res,
                                   double m1, double m2,
                                   int channelIndex = -1) const;
    //NKAON
    double NKaonScatteringXSec(double sqrt_s, const ResonanceParams& res, double massProton, double massPion) const;
    double pKaonMinus(const double &s);
    // double nKaonMinus(const double &s);

    //NPION
    // double calculateResonanceSigma(double sqrt_s, const ResonanceParams& res, double massNeutron, double massPion) const;
    // double calculateResonanceSigmaForProtonPion(double sqrt_s, const ResonanceParams& res, double massProton, double massPion, int channelIndex) const;
    // double calculateResonanceSigmaForNeutronPion(double sqrt_s, const ResonanceParams& res, double m1, double m2,int channelIndex) const;
    
    // double pPiMinus(const double & s);
    // double pPiPlus(const double & s);
    // double pPiZero(const double &s);
    // double pPiZeroWignerEckart(const double &s);
    
    // //npi
    // double nPiMinus(const double & s);
    // double nPiPlus(const double & s);
    // double nPiZero(const double & s);

    //Wrapper functions
    double XStot(const ATOOLS::Flavour & A, const ATOOLS::Flavour & B, const double & s) override;
    double XSel(const ATOOLS::Flavour & A,const ATOOLS::Flavour & B, const double & s);
  };
}








#endif