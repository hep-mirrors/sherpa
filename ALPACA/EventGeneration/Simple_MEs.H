#ifndef ALPACA_Main_Simple_MEs_H
#define ALPACA_Main_Simple_MEs_H

#include "ATOOLS/Phys/Flavour.H"

namespace ALPACA {
  class Simple_MEs {
  private:
    std::vector<double> mu2_dyn;
    
    double qiqj_qiqj(const double & hats,const double & hatt,const double & hatu);
    double qiqi_qiqi(const double & hats,const double & hatt,const double & hatu);
    double qiqbi_qiqbi(const double & hats,const double & hatt,const double & hatu);
    double qiqbi_qjqbj(const double & hats,const double & hatt,const double & hatu);
    double qqb_gg(const double & hats,const double & hatt,const double & hatu);
    double gq_gq(const double & hats,const double & hatt,const double & hatu);
    double gg_gg(const double & hats,const double & hatt,const double & hatu);
    double gg_qqb(const double & hats,const double & hatt,const double & hatu);

    double qiqj_qiqj_mu2(const double & hats, const double & hatt);
    double qiqi_qiqi_mu2(const double & hats, const double & hatt);
    double qiqbi_qiqbi_mu2(const double & hats, const double & hatt);
    double qiqbi_qjqbj_mu2(const double & hats, const double & hatt);
    double qqb_gg_mu2(const double & hats, const double & hatt);
    double gq_gq_mu2(const double & hats, const double & hatt);
    double gg_gg_mu2(const double & hats, const double & hatt);
    double gg_qqb_mu2(const double & hats, const double & hatt);

  public:
    Simple_MEs();
    ~Simple_MEs();
    
    double operator()(const ATOOLS::Flavour & in1,const ATOOLS::Flavour & in2,
		      const ATOOLS::Flavour & out1,const ATOOLS::Flavour & out2,
		      const double & hats,const double & hatt,const double & hatu);
    double operator()(const ATOOLS::Flavour & in1,const ATOOLS::Flavour & in2,
		      const double & hats,const double & hatt,const double & hatu);
    double operator()(const ATOOLS::Flavour & in1,const ATOOLS::Flavour & in2,
          const ATOOLS::Flavour & out1,const ATOOLS::Flavour & out2,
          const double & hats);
    void SetMu2Dyn(const std::vector<double> mu2_temp) { mu2_dyn = mu2_temp; }

  };
}

#endif
