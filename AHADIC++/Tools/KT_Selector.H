#ifndef AHADIC_Tools_KT_Selector_H
#define AHADIC_Tools_KT_Selector_H

#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Org/Message.H"


namespace AHADIC {
  class KT_Selector {
  private:
    //bool   m_isgluon;
    //double m_sigma,m_sigma2,m_sig2,m_offset;
    std::vector<double> m_sigma;
    std::vector<double> variation_weights;
    std::vector<double> tmp_variation_weights;
    double WeightFunction(const double & kt);
    double Gaussian(const double, const double);
    double Inv_Gaussian(const double);
    double Erf(const double, const double);
  public:
    KT_Selector();
    ~KT_Selector();
    void Init(const bool & isgluon);
    void reset_variationweights() {
      std::fill(variation_weights.begin(), variation_weights.end(), 1);
      std::fill(tmp_variation_weights.begin(), tmp_variation_weights.end(), 1);
    };

    std::vector<double>& get_variationweights() {
      return variation_weights;
    };

    void accepted() {
      for(int i{0}; i<variation_weights.size(); ++i) {
	//msg_Out() << tmp_variation_weights[i] << std::endl;
	variation_weights[i] *= tmp_variation_weights[i];
      }
    }

    double operator()(const double & ktmax,const double & M2=1.);
  };
}

#endif
